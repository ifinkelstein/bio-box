# C-21 design: a contig can belong to several genomes

Card C-21 (design half, tier H). DEC-26 (2026-11-11) let the 50 Gbp build continue by giving a shared contig to the
stage that stored it first. This note separates the two questions DEC-26 fuses, says what each option does to the
keys, joins and counts, and recommends. Nothing was run for it: it reads the code on `plan/task-cards` at
`ad4afca` and the records of measurements already made. Labels: *measured* (with its record), *derived* (formula
and inputs), *estimate* (with its basis).

**Recommendation in short.** Membership becomes many-to-many now. The rows already exist per stage in `genome_input`,
so a release reads them through a view and needs no new fact table. Gene calls stay one set per contig, and from
gene-call policy p3 they are made without genome context: a pure function of the contig's bytes and its call
inputs. p3 is adopted only if one cheap measurement (C-21f) shows the calls lose nothing that matters. Releases
built under DEC-26 are not rebuilt.

## 1. How common shared contigs become

What is shared is a whole contig, byte for byte (`contig_uid` = BLAKE3 of the sequence, `src/biobox/ids.py:10-12`).
A prophage integrated in two chromosomes does not make a shared contig. Only whole replicons do: plasmids,
phages that live as extrachromosomal replicons, and identical chromosomes. Near-identity is a different relation
(`contig_member`, spec invariant 3). Two isolates whose every contig is identical are one genome, because
`genome_uid` is the hash of the sorted contig uids (`ids.py:20-22`). Membership becomes many-to-many only when a
contig lies in two genomes with different contig sets.

| Set | Shared contigs | Label and basis |
|---|---|---|
| 50 Gbp set, one genome per species (S1, P1 to P6) | 4 of 54,265: two plasmids (P1/P3, P2/P4), two S1 phages that are contigs of P5 and P6 genomes | *measured*, C-20 probes, wilkcomp02, 2026-11-11 (DEC-26; TASKS C-20). 54,265 = 27,258 + 539 + 5,018 + 21,454 − 4 (*derived*: stage contig counts in `2026-11-10_import_ladder/README.md:17-19` and `2026-11-11_r50_inputs.md:196`) |
| INPHARED, every isolate kept | 6,818 of 34,076 records (20.0%) repeat another record's sequence | *measured*, `2026-11-10_rebuild_10gbp.md` (branch `build/r10`), inputs table. These are single-contig genomes, so they collapse to one `genome_uid` and add accessions (`contig_label.n_records`), not memberships |
| Strain collection (S2C, ≤20 complete genomes per species) | between about 4 and 3 × 10⁴, most likely of order 10³ | *estimate*. Upper bound: 56,387 contigs − 25,458 genomes = 30,929 contigs that are not a genome's first replicon (*derived*, counts *measured* in `STAGES.md:249`). Lower bound: the one-per-species rate (2 plasmid pairs in 11,966 prokaryote genomes) scaled to 25,458 genomes, about 4 (*derived*). The point estimate assumes small plasmids repeat among isolates of one species about as often as whole phage genomes repeat in INPHARED (20%), applied to an assumed half of the 30,929 being small plasmids |
| Plasmid database (S2, PLSDB, 72,360) against S2C | order 10⁴ | *estimate*. PLSDB is built from NCBI plasmid records, and the plasmid replicons of complete assemblies are among them. `STAGES.md:251` already expects "plasmid replicons of these assemblies that PLSDB already supplied at S2" to dedup by sequence uid |
| Metagenome assemblies (MAGs at S5/S6, SPIRE and Logan at S7) | rare for long contigs, unknown for short ones | *estimate*, with no number. Independent assemblies end contigs at sample-specific places, so whole-contig identity needs a complete small circle (plasmid, phage) or a short contig |

In the plasmid case, under DEC-26 and plan order (S2 before S2C), every one of those S2C genomes is one or more
contigs short in `contig_dim`. That alone means C-21 must land before S2 and S2C share a release.

What settles it: card C-21e. Pilot: the 10 GTDB R226 species with the most `Complete Genome` assemblies, ≤20 each
by the S2C rule (about 200 genomes, about 1 Gbp, *estimate*), fetched through the `ncbi_datasets_genome` form in
`sources.yaml`. Hash every contig and count the uids found in two or more genomes, by length class. Cost: about
0.3 GB of download (*derived*: 1 Gbp × 0.29 downloaded bytes per bp, `2026-11-11_r50_inputs.md:148`) and minutes
of one CPU. The full S2C count comes at no extra cost when S2C is fetched: about 31 GB in about 82 min (*derived*:
106.4 Gbp × 0.29, at the 6.28 MB/s *measured* there). PLSDB ∩ S2C can be counted by accession from metadata, but
PLSDB's URL is not in `sources.yaml` (AGENTS §8), so that needs a PI entry first.

## 2. The two questions DEC-26 fuses

### 2a. Membership: which genomes contain this contig

The facts per stage already exist. `genome_input` (`schemas.py:704-721`) holds one row per (contig, stage) with
its `genome_uid`. Ingest stores a repeated sequence once and still writes the later stage's occurrence and label
rows (`l0seq/store.py:30-34`, `567-582`). C-20 found P3's `inputs` step published before `genecall` failed, so a
later stage's `genome_input` row for a shared contig is on disk. Two things lose the second membership:

- **Within one stage** a contig in two genomes is refused, `inputs.py:57-59`. The build script hides this by
  dropping the repeat from every later genome in accession order (`scripts/stress/build_stage_pod/make_inputs.py:45-52`).
  That changes the later genome's uid and its single-mode training set, silently.
- **Across stages** `dim` refuses differing rows (`dim.py:93-97`). DEC-26 replaces the refusal with "keep the first
  stage's genome". `contig_dim` (`schemas.py:763-785`) has key (`contig_uid`, `release`) and room for one `genome_uid`.

| Option | Keys | Joins and counts | Build order | AGENTS §3.10 |
|---|---|---|---|---|
| M0, today (DEC-26) | `contig_dim` one genome per contig | no fan-out; per-genome counts miss contigs | first stage wins | kept |
| M1, `contig_dim` keyed (contig, genome) | contig- and genome-level attributes repeated on every row | every join by `contig_uid` fans out, including the ones that must not (§3) | free (set union) | kept |
| **M2, membership as its own relation** | `genome_input` key (`contig_uid`, `genome_uid`, `stage`); release view `genome_member` over the release's stages; `contig_dim` stays one row per contig; new per-release `genome_dim` | only queries that ask for genomes join it, and contig- or locus-level counts never do | free (set union) | kept: no stage file changes; new rows only in new stages and new `release=` partitions |

`genome_dim` does not exist. AGENTS §3.6 names it, but `schemas.py` has no such table. Genome-level attributes now
sit per contig in `genome_input` (`kind`, `genetic_code`, `assembly_quality`, checked equal across a genome's
contigs in `inputs.py:88-97`) and in `contig_dim` (`tier`, `assembly_quality`, `tax_order`, `is_species_rep`).
Under M2:

- `genome_dim` gets one row per genome and release: `genome_uid`, `n_contigs`, `genome_len`, `kind`, `genetic_code`,
  `assembly_quality`, `tier`, `tax_order`, `is_species_rep`.
- `contig_dim` keeps the contig-level attributes (`contig_len`, `is_circular`).
- Its `genome_uid` becomes the least uid among the contig's genomes, a function of the set.
- It gets a nullable `n_genomes` (LY-09 append; `catalog.compatible_columns`, `catalog.py:137-149`).
- Its genome-level columns are filled only where all the contig's genomes agree; otherwise they are null. A
  plasmid in hosts of two orders has no single `tax_order`, and that is the dispersion signal.

`genome_member` is a projection of `genome_input`, which is sorted by `contig_uid`. It has no aggregate, because
DuckDB pushes no join filter through a GROUP BY (`catalog.py:751-753`, *measured*). A lookup by contig therefore
reads one row group per stage file, as `contig_label` does today (*measured* 3 of 3 on `rN.1`,
`2026-11-10_rebuild_10gbp.md`). A per-release copy sorted by contig is warranted only if `profile_query` shows
this misses an SLO (spec invariant 9; card B-27 is the precedent).

### 2b. Annotation: which gene calls this contig has

Calls depend on the genome only through `_finder` (`l1/genecall.py:114-123`). A prokaryote genome of ≥100 kb is
trained on all of its contigs together. Viruses (pyrodigal-gv meta), plasmids, and genomes under 100 kb use meta
mode, and meta mode scores each sequence on its own. ncRNA hits are already independent of the other contigs
(`l1/ncrna.py:5-6`), and domains are keyed by protein. So the measured disagreements all come from single-mode
training: the *Roseovarius* plasmid has 8 genes in one host and 6 in the other, and c-st has 197 genes in gv-meta
and 198 in single mode, 11 of them different, 5.6% (*derived*, 11/198).

| Option | `feature_uid`, tapes | Counts | Build order |
|---|---|---|---|
| A1, first stage owns the calls (DEC-26) | unchanged | a locus counts once | depends on order; within a stage, undefined (refused today) |
| A2, one call set per (contig, genome context) | `blake3(contig_uid, start, end, strand, kind)` is the same where two contexts agree. `feature` (key `feature_uid`, `producer`), `contig_tapes` (key `contig_uid`, `tape_version`), `tape_tokens` and `token_postings` all need a context column. Ordinals differ per context | a locus is counted once per context, unless every query picks a context. A gene "absent" from one host is a caller artefact. Rejected | free, at the cost above |
| **A3, context-free calls (p3)** | unchanged | a locus counts once | free in content wherever the call inputs agree (§4) |

**A3 in detail.** "Call shared contigs in meta mode" cannot be decided when the first stage stores the contig, because
it is not yet shared, and changing it later would rewrite that stage. The mode must therefore depend only on the
contig. Under p3:

- viral contigs: gv-meta, as now;
- other contigs of ≥100 kb: single mode trained on the contig alone, with its genetic code;
- shorter contigs: meta.

A chromosome's training set then loses only the other replicons, and DEC-13's untestable "retrain plasmids with
the host at S2C" goes away. The c-st case shows what is at stake: gv-meta against single mode differ on 5.6% of
genes. The plasmid case shows host training moving 2 of 8 calls; my reading, not measured, is that single-mode
calls on a small foreign replicon are an artefact of the chromosome's training statistics, not information.

What to measure (card C-21f, on P1's 250 genomes and 539 contigs, read-only):

- the share of chromosome genes identical between genome-trained and self-trained calls (expected above 99%, an
  *estimate*: the training set changes by the plasmids' share of bases);
- for the non-chromosome contigs, genes matched by stop codon against the RefSeq CDS, under meta and under
  host-trained single mode, using the method behind the 95.7% K-12 figure (TASKS §1). This needs an annotation URL
  in `sources.yaml`;
- the share of bases in contigs under 100 kb.

Cost: about 1 CPU-hour (*derived*: two re-calls of P1 at the 15.39 CPU-min *measured* for its gene calling,
`2026-11-10_rebuild_10gbp.md`). Meta mode costs about 5.8 times single mode per base (*derived*: 142.79 CPU-min /
1.664 Gbp for S1's gv-meta against 143.29 / 9.596 Gbp for P2's single mode, same record). The bases share decides
whether this matters.

## 3. Every table and reader under M2 + A3 (M2 + A1 differs only where noted)

| Table | Schema, key, sort key, lock | What changes |
|---|---|---|
| `genome_input` (`schemas.py:704-721`) | key → (`contig_uid`, `genome_uid`, `stage`), sort → (`contig_uid`, `genome_uid`); columns and lock unchanged; older files still sorted | a contig is listed in every genome of the stage that holds it. Disagreeing `kind`, `genetic_code` or `is_circular` within one stage is refused, naming both genomes |
| `contig_occurrence`, `contig_label`, `contig_index` (`:90-127`, `:359-380`) | none | already record-level, or written once by the first stage |
| `contig_dim` (`:763-785`) | append nullable `n_genomes`; lock regenerated | `genome_uid` = least member; genome attributes only where the genomes agree |
| `genome_dim` (new, `dim/genome_dim`, `release`) | key (`genome_uid`, `release`), sort `genome_uid`; lock regenerated | written by `dim` |
| `feature`, `feature_protein`, `protein_occurrence`, `igr`, `contig_tapes` (`:722-762`, `:217-239`, `:159-179`, `:630-665`) | none | written by the first stage that stores the contig. A later stage writes nothing for it (C-20). Under A3 its own call would be byte-identical (a test, not a write); under A1 it would differ |
| `tape_tokens`, `token_postings` (`:683-702`, `:821-850`) | none | one tape per contig, so one posting per slot and token, so a locus counts once however many genomes carry it |
| `l1_domains` (`:240-272`) | none | per protein |

How each reader handles a hit on a shared contig. Counts by genome fan out; counts by contig and locus must not.

- **`q1_count`** (`catalog.py:895-899`) counts distinct (contig, ordinal); keep that. Its `tier` filter moves to a
  semi-join, `contig_uid IN (SELECT contig_uid FROM genome_member JOIN genome_dim_now() ... WHERE tier = _tier)`:
  a semi-join cannot fan out. "Number of genomes with this neighbourhood" is a new macro, `q1_genomes(a, b, k)`,
  that counts distinct `genome_uid` over the loci's contigs joined to `genome_member`.
- **`pairs_of`, `fusions_of`, `tapes_of`** are per contig and unchanged.
- **`partners_of`** (`catalog.py:976-995`) keeps `n_loci` and `n_pairs` as locus counts. A partner's genome count is
  the union of the genomes of its loci's contigs. It can be an appended `n_genomes` column, computed in
  `bb_partner_counts` from a pruned `genome_member` lookup of the anchor contigs, which is the same bound as its
  `tape_tokens` read.
- **New macros:** `genome_dim_now()`, and `genomes_of(_cuid)`, a pruned lookup.
- **`search hits` (hits.tsv)** keeps one row per gene (`export.py:555-575`). `genome_uid` is `contig_dim`'s, which is
  already one per contig (`export.py:500-501`). It adds an `n_genomes` column. The JSON `search_hit` occurrences
  (`export.py:88-104`) add `n_genomes` and `genome_uids` (key additions need no version bump, `export.py:237-252`).
  Each loci file gets `n_genomes=` and `genome_uids=` in its facts line. The FASTA is unchanged.
- **`nbr summarize`** breaks as soon as `genome_of` has two rows for one contig. `_aggregate` (`nbr.py:357-377`)
  joins `genome_of` to `nbr_rows` and to `_loci` with LEFT JOIN, then counts `count(*)`. That would double the
  per-key `n_neighbors` and offset counts, and the header's `n_neighbors`, `n_anchor_loci` and `n_censored_loci`.
  Per-key `n_loci` is a distinct count and survives, and so does `same_strand_share`, a ratio of two doubled counts.
  The fix:
  `genome_of` reads `genome_member` (`nbr.py:533`), and genome counts are computed as `count(DISTINCT genome_uid)`
  in a separate aggregate that never joins into the row-level tables. `n_loci`, `n_contigs` and `n_neighbors` stay
  by locus.

## 4. Order independence

- **M2** is a pure function of the release's set of stages: a union of `genome_input` rows, `least()` and `count`.
- **A1** gives a release that is a function of its stages in build order. Within a stage there is no order.
- **A3** gives content that does not depend on order wherever the genomes agree on a contig's call inputs (caller
  class, code, topology). What still depends on order is the layout: which stage's partition holds the rows. That
  is exactly what DEC-17/21 already accepted for sequences (`seq_locator` per stage, first stored), and the
  contig's `contig_index` row records it.
- **The residue under A3** is contigs whose genomes disagree on call inputs. Two of the four measured are of this
  kind: the S1 phages, which are `virus` in S1 and `prokaryote` in P5 and P6. They keep the first stage's inputs
  and are written to the ledger as `call_inputs_conflict`, with a count in the release report.

**No option is order-free everywhere without rewriting earlier stages.** Even A3 is order-free only once call
inputs are functions of the bytes, which needs a contig-level classifier (geNomad, already in spec invariants 3
and 8) and a topology rule. The cheapest honest statement for the release report:

> Membership depends only on the stages. A shared contig's gene calls are a function of its sequence and call inputs (p3), and of the earliest stage in build order where those inputs disagree (N contigs, listed).

Under A1 the statement instead reads "the calls of the earliest stage that contains it, or of the genome with the
least uid within a stage". It is still honest, but it covers every shared contig, which could be 10³ to 10⁴ by S2C
(estimate, §1).

Is A3 worth it? It costs no migration of production data: `rS0.1` is not signed, and the 10 and 50 Gbp roots are
pilots. It adds meta-mode CPU on short contigs, and it changes calls on plasmids, which C-21f measures. That buys
order-free calls for strain collections and for the S2/S2C overlap, and removes the within-stage owner choice. I
think it is worth it if C-21f shows that chromosome calls agree to within about 1% and plasmid calls are no worse
against the reference. Otherwise, keep A1 with the statement above.

## 5. Migration

- **`rN.1` to `rN.3`** (stages S1, P1, P2) contain no shared contig, because each measured pair has a member in
  P3 to P6 (*derived* from the C-20 pair list). They are answered as today.
- **The 50 Gbp release** holds the 4 shared contigs.
  - Its calls stay DEC-26's.
  - Its membership is answered without a rebuild: `genome_member` is a view over `genome_input` files that the
    release's manifest already lists. For a catalog built before the change, `open_catalog` creates the view and
    macros as temporary objects, as it already does for `partners_of` (`catalog.py:929-933`).
  - That release's `contig_dim` lacks `n_genomes` (null), so readers count from the view.
  - Check: C-20's code must keep the later stage's `genome_input` row (open item 5).
- **Nothing is rewritten.** New stages write the new `genome_input` key, and new releases write `genome_dim` and
  `n_genomes`.
- **p3** is a new `GENE_CALL_VERSION` (`genecall.py:77`). It cannot be added beside p2 in an existing root:
  `l1_domains` is partitioned by (`library`, `stage`) without the gene-call version (`schemas.py:263-265`), and a
  release serves one tape version. So p3 starts in a new root, which should be the production root for S0 onward.

## 6. Recommendation, decision rows, cards

Recommendation: M2 now, before S2C. A3 (p3) before the production S0 build, if C-21f passes; otherwise A1 with
the release-report statement of §4.

Draft rows (numbers to be assigned by the PI):

> **DEC-27** | **A contig may belong to several genomes (amends DEC-26).** `genome_input` lists a contig once for each genome that contains it, also within a stage, key (`contig_uid`, `genome_uid`, `stage`). A release's membership is the union of its stages' `genome_input` rows, read through the catalog view `genome_member`. `contig_dim` keeps one row per contig, with `genome_uid` the least uid of its genomes and an appended nullable `n_genomes`. Its genome attributes are null where its genomes disagree. Genome attributes live per release in `genome_dim`. A count by genome counts distinct genomes through `genome_member`; a count by contig, locus or neighbour never joins it except as a semi-join. Within a stage, a shared contig's calls are written by the genome with the least uid. Releases built before this keep their files. Basis: C-20 (4 of 54,265 measured; estimate 10³ to 10⁴ at S2/S2C, C-21 note §1). | `planning/design/C-21_shared_contigs.md`; `src/biobox/l0seq/inputs.py`, `src/biobox/dim.py`, `src/biobox/catalog.py`, `src/biobox/search/` |
>
> **DEC-28** (after C-21f) | **Gene calls belong to the contig, not the genome: policy p3 (amends DEC-13, DEC-26).** Viral contigs: pyrodigal-gv meta. Other contigs ≥100 kb: single mode trained on the contig alone with the genome's code. Shorter contigs: meta. A shared contig's calls are written once, by the first stage that stores it, and are the same whichever stage that is. Genomes that disagree on a contig's caller class, code or topology keep the first stage's inputs, recorded as `call_inputs_conflict`. Starts in a new root, with a new `gene_call_version`. Basis: C-21f. | `src/biobox/l1/genecall.py` |

Cards, each **M**, test first (AGENTS §5.1):

- **C-21a `genome_input` lists a shared contig in every genome** (depends on C-20).
  - Files: `src/biobox/l0seq/inputs.py`, `src/biobox/registry/schemas.py` (`genome_input` keys only),
    `src/biobox/l1/genecall.py`, `src/biobox/l1/ncrna.py`, `src/biobox/stage.py`,
    `scripts/stress/build_stage_pod/make_inputs.py`, tests.
  - Failing test: a stage of two genomes sharing a plasmid raises at `inputs.py:59`.
  - Then: both rows are written; the plasmid's calls are written once, by the genome with the least uid; the output
    is byte-identical with the genomes given in either order; disagreeing call inputs within a stage are refused;
    `make_inputs` stops dropping repeats and reports them.
  - Accept: the tests; the lock unchanged.
- **C-21b membership view and `genome_dim`** (depends on C-21a).
  - Files: `src/biobox/dim.py`, `src/biobox/registry/schemas.py` (`genome_dim`; `contig_dim.n_genomes`, regenerated
    lock), `src/biobox/catalog.py` (`genome_member`, `genome_dim_now`, `genomes_of`, and their creation in
    `open_catalog` for old catalogs), tests.
  - Failing test: after C-20, a release over two stages that share a contig reports it in one genome only; the
    test expects two rows in `genome_member` and `n_genomes` 2.
  - Then: `contig_dim`, `genome_dim` and `genome_member` are byte-identical when `release_def` lists the stages in the
    other order; `genomes_of` with 100 joined keys reads at most one row group per key and stage file on the 100×
    synthetic store (`profile_query`); a memory test at the card's bound.
- **C-21c counts by genome fan out, counts by locus do not** (depends on C-21b).
  - Files: `src/biobox/catalog.py` (`q1_count` tier semi-join, `q1_genomes`, `partners_of` `n_genomes`),
    `src/biobox/search/nbr.py`, `src/biobox/search/export.py` (README text included), the W-07 byte-identity test
    and the JSON schema lock (regenerated).
  - Failing test, a metamorphic relation: add a stage whose genome repeats an existing contig. Today the genome
    counts on that contig stay at 1 (`contig_dim` has one genome). Pointing `genome_of` at `genome_member` without
    the fix doubles `n_neighbors` (§3), which is the second failing assertion.
  - Then: postings, `q1_count`, `pairs_of`, `partners_of`'s `n_loci` and `n_pairs`, `nbr` `n_loci`, `n_neighbors`,
    `offsets` and the hits.tsv row count are unchanged, and every `n_genomes` on that contig rises by one.
- **C-21d gene-call policy p3** (depends on C-21a and DEC-28).
  - Files: `src/biobox/l1/genecall.py`, `src/biobox/stage.py`, tests.
  - Failing test: the same plasmid in two genomes of different composition gives different features.
  - Then: identical rows; per-contig calls equal the Prodigal 2.6.3 binary (`-p single` on the contig alone, `-p
    meta`) as the external oracle; conflicts are written to the ledger.
- **C-21e measure shared contigs in a strain collection** (Pod: yes; §1 pilot) and **C-21f measure p3 call
  quality** (Pod: yes; §2b). Each writes `planning/measurements/<date>_*.md` and changes no `src/` code. C-21f needs
  a `sources.yaml` entry for reference annotations first.

## 7. Not determined here, and the measurement that settles each

1. **Shared contigs in a strain collection, by length class.** Settled by C-21e, the 200-genome pilot.
2. **PLSDB ∩ S2C.** Settled by an accession join of the PLSDB and S2C metadata (needs PLSDB in `sources.yaml`).
3. **p3 quality: chromosome agreement, and plasmid calls against reference CDS.** Settled by C-21f.
4. **Share of prokaryote bases in contigs under 100 kb** (the meta-mode CPU of p3). Settled by one query over the
   50 Gbp root's `contig_index` and `genome_input`, seconds.
5. **Whether C-20's code keeps the later stage's `genome_input` row for a shared contig.** Settled by one count on
   the 50 Gbp root after P3 completes: rows whose `contig_uid` has no `contig_index` row in their own stage. 4 are
   expected.
6. **Exact sharing among MAGs.** Settled by hashing the genomes of one GTDB species with many MAGs.
7. **How often genomes disagree on call inputs** (the residue of §4). Settled by items 1 and 2 above, by
   comparing `kind`, `genetic_code` and topology across a shared contig's genomes.
8. **Row groups per key for `genome_member` at S6's stage count.** Settled by `profile_query` on the 100× synthetic
   store split into 20 stages.
