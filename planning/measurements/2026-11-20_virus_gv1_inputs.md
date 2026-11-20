# GVMAGs V2 as stage GV1: inputs, overlap and cost (card C-48c, part 3)

Date: 2026-11-20, wilkcomp02, 14:23 to 14:42 CST (the pod's clock). Every number carries a label: *measured* (host,
date, command or file; wilkcomp02, 2026-11-20 unless said), *derived* (formula and inputs), *estimate*, or *published*.
**No build**: nothing was written in the data root `/tmp/if2453/r10/root` (bound read-only) or under `/stor` (the
archive bound read-only); V1's and V2's inputs bound read-only. No external fetch.

Part 2 (`2026-11-20_virus_v2_gv1_inputs.md`) measured GVMAGs V2 and stopped GV1 at its selection rule. The
orchestrator's decisions for this part, each about inputs only: (1) GV1 is **every** genome of
`all_18kgenomes_fna.tar.gz` (18,726), since the archive has no quality table to select by and DEC-59 takes everything in
the end; (2) V2 is built before GV1, so V2 owns the contigs the two share and GV1 lists them as second genomes
(DEC-26, DEC-27); (3) assembly quality is the least committal value the registry allows; (4) no terminal-repeat
trimming, every contig as the archive gives it; (5) kind virus, gene-call mode gv-meta.

Beside this file in [`2026-11-20_virus_gv1_inputs/`](2026-11-20_virus_gv1_inputs/): `inputs_SUMS`, `report.json`
(staging), `check.json` (dry check), `overlap.json`, `gvsubset_summary.json` (fna_files against GV1), `est.json`,
`objects.json`, `repeated_records.tsv`, `shared_records.tsv`, `skipped_members.tsv`, `run.log`. Scripts:
`c48c_virus_inputs.py gv1` (tests `tests/r1/test_c48c_virus_inputs.py`), `c48c_v2_gv1_overlap.py --gv1-inputs`,
`c48c_v2_gv1_look.py gvsubset`, `c48c_v2_gv1_est.py stage`, run by `c48c_v2_gv1_run.sh gv1-<mode>` in units
`c48c3-gv1-<mode>` (pod tree `/tmp/if2453/c48c3/tree`, log `/tmp/if2453/c48c3/logs/run.log`, outputs
`/tmp/if2453/c48c3/work/<mode>/`).

## Summary

- **Stage GV1** is staged at `/tmp/if2453/s2/inputs/GV1/`: **18,499 genomes from 18,726 `.fna` members** (227 repeat
  an earlier genome's contig set exactly, DEC-27), **281,804 distinct contigs, 4,316,061,173 bp**, kind virus
  (gv-meta), code 11, quality `contig`, every contig linear, CC BY-SA 4.0 (redistributable) (*measured*,
  `report.json`). The archive member's md5 equals `sources.yaml`. Dry check exit 0, no problem; a second run gives
  every file byte for byte.
- **Shared**: 17,670 contigs (425,418,119 bp) are held by an earlier stage: **V2 17,664** (its IMG/VR GVMAG
  representatives) and **H1 6**; none by S1, P1 to P6, PL1, C5 or V1. **GV1 adds 264,134 contigs, 3,890,643,054 bp.**
  1,116 of GV1's genomes have every contig held earlier: they stay genomes of GV1, and the build calls none of them.
- **Cost** (*derived*): `annotated` 9,471 CPU-min at the code's rates (Pfam 9,076 in 38 shards); `called` 382. Gene
  calling measured here at 2.72 s per Mbp, **35% under** the code's gv-meta rate: at the measured rate the gene-calling
  step is 230 CPU-min, not 354, so `called` is 258 and `annotated` 9,347. Proposal below (§5).
- **Nothing of GVMAGs V2 is left for later**: 18,726 members = 18,499 staged + 227 repeated; the one other member is
  an editor backup of a shell script. The `fna_files` subset (8,508) is the same genomes by name: 8,066 with the same
  contigs and 442 with 746 contigs (5.3 Mbp) fewer; it holds **no contig GV1 lacks**.

## 1. The rules (test first)

`c48c_virus_inputs.py gv1 GVMAGS_TAR OUT --root ROOT --root-stages S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5 --owner-inputs
V1=/tmp/if2453/s2/inputs/V1 --owner-inputs V2=/tmp/if2453/s2/inputs/V2`:

- The member `gvmags_v2/all_18kgenomes_fna.tar.gz` is streamed out of the outer tar through a hashing reader (its md5
  checked against `sources.yaml`'s Zenodo md5, `077f125d…40b9`) and read one `.fna` member (one genome, a few MB) at a
  time; only per-record uid, length and name are held (AGENTS §4.1).
- **A genome is one `.fna`**: its records are its contigs, in file order. Kind virus, so gv-meta mode
  (`genecall.mode_of`; pyrodigal-gv's giant-virus models). Genetic code 11 recorded (gv-meta takes each gene's code
  from its models).
- **Assembly quality `contig`**: of the schema's values (complete, chromosome, scaffold, contig, mag_hq, mag_mq, vOTU,
  plasmid), `mag_hq` and `mag_mq` claim a completeness the archive does not state, `complete`, `chromosome` and
  `scaffold` an assembly level, `vOTU` a clustering; `contig` claims nothing beyond "assembled contigs", and it is
  the value card C-40 gives a plasmid of unknown completeness (`make_inputs.COMPLETENESS["unknown"]`). No change under
  `src/` was needed.
- **Every contig linear**: the archive records no topology; no repeat is trimmed (decision 4).
- A genome whose contig set an earlier genome has is that genome (listed once; `repeated_records.tsv`). A contig in
  several genomes of different contig sets stays in each (DEC-27): `genome_input` lists it once per genome, and
  gene calling writes its rows once, from the genome with the least `genome_uid` (`genecall`'s DEC-27 rule).
- A member that is not `.fna` is not read as a genome and is listed in `skipped_members.tsv` (never run).
- Owners as for V2: the root's stored stages' `contig_index` default variants, then V1's inputs, then V2's inputs;
  a record whose contig one holds stays in its genome and is written to `shared_records.tsv` with the first owner.

Tests (`tests/r1/test_c48c_virus_inputs.py`, two new; synthetic tar in `tmp_path`): a genome of three contigs one of
which V2 holds; a genome both of whose contigs V2 holds; a genome with another's bytes under another name; a genome
sharing one contig with another; a non-FASTA member; the md5 refusal; two runs byte for byte. They failed first on the
pod (`2 failed, 9 passed`: `AttributeError: module 'c48c_virus_inputs' has no attribute 'stage_gv1'`); now the file's
11 tests and `tests/r0` pass (235 passed, *measured*). The part 2 tests pass unchanged.

## 2. Stage GV1's inputs

Unit `c48c3-gv1-inputs`, 14:23:06 to 14:26:13, one thread (*measured*). Written as `GV1.tmp`, renamed when complete.

| File (`/tmp/if2453/s2/inputs/GV1/`) | Content | bytes | sha256 |
|---|---|---|---|
| `genomes.json` | 18,499 genomes, 287,159 memberships, kind virus, code 11, quality contig, none circular | 12,224,622 | `a86ce75b88da4c31e7016b0c90e6b6967a1d525604d8c5247d94bfeeec595ced` |
| `objects.json` | one object, source `gvmags_v2`, zenodo:10.5281/zenodo.20126600, CC BY-SA 4.0, redistributable, retrieved 2026-11-19T10:42:14Z (the end of its fetch, `sources.yaml`) | 310 | `2eaadb35ecc8eec0d12fcbcb69cccdb2ae48758506f53383426caa3e066b18fb` |
| `records.tsv` | every record: ordinal, genome, record, accession, contig uid, length, member | 39,879,303 | `a87ceb6981f9cce20994eba471a01088a31ba5cff1b8f497c7d677ccdf321380` |
| `repeated_records.tsv` | 227 genomes whose contig set an earlier genome has, with that genome | 17,779 | `16f2902dcd2d3cc197abe7adb20c783a86808739b464a49de3dd1ec6c5eb32e7` |
| `shared_records.tsv` | 21,984 records whose contig an earlier stage holds, with the owner | 1,450,996 | `4fb544d58a1307071c51a6fc595e78273382fc8c36b269b842d7557db72cdb88` |
| `skipped_members.tsv` | `confirmed/chunk_maker.sh~` (1,001 bytes), not_fasta | 61 | `2679a24e938245a54dedd841d7bbf4369bcb4e3cd22a8f19620a393e7cabfad3` |
| `gvmags_v2_all_18kgenomes.fna` | 288,536 records, 4,576,307,473 bp, headers and line width as in the archive, archive order | 4,639,876,745 | `568f3a312fd78d03544792bfb3a89bde708591865065d72c62fc6c0dbdc0028b` |

*Measured* (`inputs_SUMS`, `report.json`).

- **Dry check** (`s2_check_inputs.py GV1`, unit `c48c3-gv1-check`): exit 0, 18,499 genomes, 287,159 `genome_input`
  rows, 18,499 distinct genome uids, 0 circular rows, kinds `[virus]`, codes `[11]`, quality `[contig]`, gene-call
  mode gv-meta, `problems: []` (*measured*, `check.json`).
- **Determinism**: unit `c48c3-gv1-rerun` (14:29:28 to 14:32:47) gave all seven files byte for byte (`cmp`;
  `objects.json` after its own path) (*measured*, `run.log`).
- **Count check** (*measured*): 18,726 `.fna` members = 18,499 genomes listed + 227 repeated; 1 non-FASTA member
  skipped; no member without a record; no contig repeated inside a genome.
- **Contigs in several genomes**: 5,255 contigs are listed by more than one of the 18,499 genomes (*measured*,
  `report.json`). Part 2's 6,626 counted genome names, the 227 repeated genomes among them. Each is one stored contig
  with one call set (the least `genome_uid` writes it); every genome holding it lists it (DEC-27). Median contigs per
  genome 11.

## 3. What GV1 shares by `contig_uid`

`c48c_v2_gv1_overlap.py --gv1-inputs` (unit `c48c3-gv1-overlap`, 14:29 to 14:29:41), DuckDB on each stored stage's
`contig_index` and `genome_input`, V1's and V2's `records.tsv` and GV1's (*measured*, `overlap.json`; the same counts as
part 2's measurement on the archive):

| | S1 | P1–P6 | H1 | PL1 | C5 | V1 inputs | V2 inputs | held earlier | **adds** |
|---|---|---|---|---|---|---|---|---|---|
| GV1 contigs | 0 | 0 | 6 (2,242,617 bp) | 0 | 0 | 0 | 17,664 (423,175,502 bp) | 17,670; 425,418,119 bp | **264,134; 3,890,643,054 bp** |

- Owners in `shared_records.tsv` (*measured*): V2 17,664 contigs in 21,977 records, H1 6 contigs in 7 records.
- 17,567 of the 17,664 contigs shared with V2 are contigs of V2's 1,422 IMG/VR GVMAG representatives (Schulz et al.
  2020's bins, part 2); the rest are 97 single-record V2 representatives (Linear 51, ITR 26, DTR 20).
- **Genomes all of whose contigs are held earlier: 1,116** of the 18,499 (*measured*, `report.json`; 1,174 counting the
  227 repeated names, `overlap.json`); 682 genome names hold some. They remain genomes of GV1 (a `genome_input` row
  per contig); `genecall` calls none of the 1,116, since every row would be dropped (its item 5), and records the
  skips.

## 4. `fna_files.tar.gz` against GV1

`c48c_v2_gv1_look.py gvsubset` (unit `c48c3-gv1-subset`, 14:29 to 14:29:58), each of fna_files' genomes' contig uids
against the GV1 genome of the same name (*measured*, `gvsubset_summary.json`):

- 8,508 genomes, every one named in GV1: **8,066 with exactly the same contigs**, **442 a strict subset** of their GV1
  genome. The 442 lack 746 of GV1's contigs (5,323,061 bp; median 7,123 bp per genome): 311 lack one contig, 64 two,
  the rest 3 to 15.
- **fna_files holds no contig GV1 lacks** (0). So the subset is the 18k set's genomes with a few contigs removed (by
  the provider, for a reason the archive does not state, perhaps decontamination); staging GV1 from the 18k set loses
  nothing of it. Whether those 746 contigs should be flagged is a question for the PI (§7).

## 5. What building GV1 would cost

`c48c_v2_gv1_est.py stage` (unit `c48c3-gv1-est`, 14:32 to 14:42:10, one thread): the code's own functions over GV1's
inputs; a sample of every 25th genome with a contig of its own (695 genomes, 174,342,637 bp) called with the step's
own path (`genecall._call`, one thread) (*measured* sample, *derived* totals, `est.json`).

| Step | at the code's rates | at the measured gene-calling rate | How |
|---|---|---|---|
| ingest | 17.4 (19.1 counting every byte) | same | `stage.ingest_cpu_min(4,639,876,745 bytes, 425,418,119 held)` |
| gene calling, gv-meta | **354.1 in 2 batches** | **229.7** (one batch would fit) | code: `one_process_cpu_min` over 4,108,278,835 bases read (every genome with a contig of its own, whole); measured: (15 s + 1.23 × 2.7239 s/Mbp × 4,108.28 Mbp) / 60 |
| genes; distinct proteins; residues | 4.13 M; 3.97 M; 1,140.3 M | same | sample per Mbp (1,060.7; 1,020.0; 293,076) × 3,890.64 Mbp stored (distinct: an upper bound) |
| Pfam 38.2 scan | 9,076.4 in 38 shards | same | 1,140.25 M residues × 7.96 (`SCAN_CPU_MIN_PER_M`) |
| CRISPR, default and near-Cas | 10.9 and 12.5 | same | `stage.crispr_cpu_min(3,890,643,054)`; no circular contigs |
| **`--depth annotated`** | **9,471.4** | **9,347.0** | the sum |
| **`--depth called`** | **382.4** | **258.0** | ingest, calls, default CRISPR |
| what `called` saves | 9,089.0 (96%) | 9,089.0 | Pfam and the near-Cas producer |
| files in the root | about 2.8 GB | | 3.89 Gbp × 0.7175 to 0.7316 GB per Gbp (card C-39) |

**The sample** (*measured*): 184,923 genes (1,060.69 per Mbp), 177,820 distinct proteins (1,019.95 per Mbp; 322 in
the root's `protein_seq`), 51,095,667 residues; calls 474.9 CPU-s, **2.7239 s per Mbp**; every sampled genome
completed. Part 2's sample (every 50th genome, 374 genomes) gave 2.6053.

**Against the code, and what a build driver should judge against.** The code's gv-meta linear rate is 4.2019 s per Mbp
(`genecall.CALL_S_PER_MBP`, measured on phages of stage S1); giant viruses run at 2.7239, **35.2% under** it
(*derived*: 229.7 / 354.1 = 0.649 on the whole step). What this changes, *derived*:

- At `annotated` the Pfam scan is 96% of the stage, and the total moves by 1.3% (9,471 to 9,347): the code's estimate
  is within AGENTS §6's 20% and a build may be judged against it.
- At `called`, and for the gene-calling step alone, the code's estimate is 48% above (382 against 258; 354 against
  230): **more than 20% off**. The code's figure is safe for planning (batches smaller than needed, never a process
  over the budget), but a driver judging the step's run record against it would accept a step 30% slower than
  expected. Proposed, for a code card (not changed here): a gv-meta rate for giant-virus genomes of **2.80 s per Mbp**
  (the largest of the two samples, 2.7239, rounded up as the code rounds its rates), which puts GV1's gene-calling
  step at **236 CPU-min** (*derived*: (15 + 1.23 × 2.80 × 4,108.28) / 60 = 236.1) and `called` at about 264; until it
  exists, a driver should judge GV1's gene-calling run record against 236, and plan batches with the code's 354.
  The Pfam rate (7.96 per M residues) has not been measured on giant-virus proteins; the first shards' run records
  should be compared with it before the rest are launched.

No single process passes 300 CPU-minutes: gene calling in 2 batches of at most 240 at the code's estimate, each Pfam
shard at most 240, ingest 19.1.

## 6. What is left for later

Nothing of GVMAGs V2 (*measured*): of `all_18kgenomes_fna.tar.gz`'s 18,727 members, 18,499 genomes are staged, 227 are
repeats of a staged genome (their records are in the FASTA object and are ingested as occurrences), and 1 is a shell
script backup (`skipped_members.tsv`): 18,499 + 227 + 1 = 18,727, of which 18,726 are genomes. `fna_files.tar.gz`
adds no contig (§4); `long_ortho_8000.tsv` (orthogroups over 8,501 genomes) is a table, not sequence, and is not
staged. Part 2's `GV1_left_for_later.tsv.gz` (18,726 genomes, reason `gv1_selection_undecided`) is superseded.

## 7. Open questions

1. **The 746 contigs** that the provider's `fna_files` subset removes from 442 genomes (5.3 Mbp): GV1 keeps them, as
   the 18k set has them. Flag them in a later dimension card, or leave them?
2. **A giant-virus gene-calling rate** (§5) for a code card: 2.80 s per Mbp proposed.
3. **Ownership of the 17,664 shared contigs** follows the build order (V2 first). If GV1 were built first, V2's 1,422
   GVMAG representatives would be the second genomes instead (part 2 §4).

## 8. Pod use and scratch

| Unit | CST | What | Threads |
|---|---|---|---|
| `c48c3-gv1-inputs` | 14:23:06 to 14:26:13 | GV1's inputs | 1 |
| `c48c3-gv1-rerun` | 14:29:28 to 14:32:47 | determinism | 1 |
| `c48c3-gv1-check`, `-gv1-overlap`, `-gv1-subset` | 14:29 to 14:30 | dry check, overlap, fna_files | 1 each |
| `c48c3-gv1-est` | 14:32 to 14:42:10 | estimates, gene-calling sample | 1 |
| `pod_pytest.sh C-48c-3` | 14:25 | `tests/r1/test_c48c_virus_inputs.py` and `tests/r0`: 235 passed | 1 |

Every unit: `systemd-run --user`, `nice 19`, `ionice` 7, `MemoryMax` 16 GB, runtime image `bio-box_faec1f4.sif`,
`PYTHONPATH` the tree's `src` alone; the data root, the archives, V1's and V2's inputs bound read-only; at most 4
threads of ours at once. No `STOP` file. (The run script logs part 3's units with its part 2 prefix `c48c2:`; the
units and paths are part 3's.) Disk (*measured*, `du -sb`, 14:43): GV1's inputs 4,693,449,816 bytes; scratch
`/tmp/if2453/c48c3` 46,511,649. Nothing on `/stor`.

## 9. Proposed register rows (for `planning/CATALOG_GOAL.md` §4.1)

| Source | Taken | Left for later, in full | Where the remainder is recorded | How to pick it up |
|---|---|---|---|---|
| GVMAGs V2 (`gvmags_v2`) | GV1: every genome of `all_18kgenomes_fna.tar.gz`, 18,499 genomes (227 repeats listed once), 281,804 contigs, 4.32 Gbp (3.89 Gbp new after V2) | none; `fna_files.tar.gz` (8,508 genomes) holds no contig GV1 lacks | `2026-11-20_virus_gv1_inputs.md` §4, §6; `2026-11-20_virus_gv1_inputs/repeated_records.tsv`, `skipped_members.tsv` | — |

## 10. Proposed `sources.yaml` lines (not applied)

Under `gvmags_v2` (replacing part 2's `staged: none`):

```yaml
    staged:
      - stage: GV1
        object: /tmp/if2453/s2/inputs/GV1/gvmags_v2_all_18kgenomes.fna
        bytes: 4639876745
        sha256: 568f3a312fd78d03544792bfb3a89bde708591865065d72c62fc6c0dbdc0028b
        derived_from: gvmags_v2_20261118.tar, member gvmags_v2/all_18kgenomes_fna.tar.gz (md5 077f125d844c85ca6169c331d52640b9),
          every .fna member's records in archive order
        rule: c48c_virus_inputs.py gv1 (every genome; kind virus, gv-meta; quality contig; card C-48c part 3)
        left_for_later: none (fna_files.tar.gz is a subset by contigs, 2026-11-20_virus_gv1_inputs.md §4)
```
