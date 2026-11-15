# Label files for the stages ingested before card B-27 (card C-34)

Stages S1, P1 and P2 of the data root `/tmp/if2453/r10/root` were ingested before card B-27 and had no
`contig_label` file; P3 to P6 have one. Since card W-36 readers take those three stages' labels from their
`contig_occurrence`. This card adds `bio-box stage labels --stage S`, the rule that lets a published release verify
after its stages gain such files, and the label files of S1, P1 and P2 in the real root.

Every number here is *measured* on wilkcomp02 on 2026-11-15, in the runtime image `bio-box_faec1f4.sif` (sha256
`d67300aa…`), with the card's tree copied to `/tmp/if2453/c34/tree`, unless it is labelled otherwise. Times are the
pod's clock (CST). The driver is `scripts/stress/build_stage_pod/c34_driver.sh`, with `c34_view.py`, `c34_oracle.py`
and `c34_check.py`; its log is `/tmp/if2453/c34/driver.log`.

## The step

`bio-box stage labels --stage S` (`biobox.l0seq.store.write_labels`) reads the stage's own `contig_occurrence` files
(every `source=` partition of `l0seq/contig_occurrence/stage=S/`).

- A DuckDB view of the label columns over those files stands for ingest's private scratch table. Ingest's code then
  writes the table: `_label_batches`, one DuckDB sort by (`contig_uid`, `l0_object_sha256`, `record_ordinal`) that
  spills under `DUCK_MEMORY_LIMIT`, read a batch at a time, and `io.write_parts` at `EST_LABEL_ROW`.
- Every record of the stage is labelled, a contig another stage owns (DEC-26) included, as ingest labels it.
- It writes the ledger `ledger/stage=S/step=labels`, one row per source object: records in; kept, the records that
  are their contig's first (the label rows naming the object); merged, the others. The tool versions are pyarrow and
  DuckDB.
- It refuses a stage that has label files or the step's ledger, in the root or in `--out`. It also refuses a stage
  whose `contig_occurrence` holds another number of records than its ingest ledger counts, and a stage with no
  record.
- Under `run-step` it is the step `labels`, with a run record.

## The rule for a release published before the files (`verify_manifest`)

**Before the card.** `contig_label` has no partition but its stage, so card C-19's `_owner_key` gave its files no
owner key. Every release listing the stage owned them, and each new label file was `not in manifest`.

On a view of the real root with the rehearsal's three label files in place (`/tmp/if2453/c34/verifyview`: links, and
real directories for `l0seq/contig_label/stage=S1|P1|P2` holding links to the files), `release verify` by the code
before the card (`60d4e7e`'s `src`, `/tmp/if2453/c34/oldtree`) failed **all six** releases with a manifest, not the
four the card names. `rN.1a` (S1) and `rN.1b` (S1, P1) are affected too:

| Release | Problems, code before the card | With the card's code |
|---|---|---|
| `rN.1a` | 1: `not in manifest: l0seq/contig_label/stage=S1/part-0.parquet` | verified, hash `39020b27de1fa163570c4255365aebdd` |
| `rN.1b` | 2 (S1, P1) | verified, `d3212c9f264e25208158534e61e24011` |
| `rN.1` | 3 (S1, P1, P2) | verified, `84c2c057b09643bb4e9925938703095d` |
| `rN.2` | 3 | verified, `2b4b0e43cecdb3a5dbccd1a2bb2b0393` |
| `rN.3` | 3 | verified, `eee598acd16324d1dd73349cfd6d931f` |
| `rM.2` | 3 | verified, `2239df7b92ec775058fb9c2ddae3609e` |

The logs are `/tmp/if2453/c34/rehearsal-verify/{old,new}-<R>.log`, from 03:05 to 03:17.

**The rule** (`catalog.STAGE_ADDED_TABLES`, used in `_owner_key`). A `contig_label` file under stage S is a release's
concern exactly when its manifest lists a `contig_label` file of S. The manifest alone tells the two cases apart:

- **It lists none for S.** The release was built before S had labels, its readers take S's labels from
  `contig_occurrence` (AGENTS §3.9: they read the manifest's lists), and every label file of S was added after it.
  Those files are ignored.
- **It lists one or more for S.** They are checked as every listed file is: `missing:`, `size differs:`,
  `sha256 differs:`. Any other file under `contig_label/stage=S/` is `not in manifest:`, as before.

A changed or missing file of any other table of S, and a path no table claims, are reported as before. No table but
`contig_label` is named: the code has no general notion of a stage table added later, beyond C-19's partitions.

**Tests.** `tests/r1/test_c34_verify_rule.py` tests both cases on a small world built in the test:

- release A lists S0's labels and not S1's;
- release B lists none;
- release C, defined after S1 gained labels, lists them.

## Tests (fixture)

The tests were written first. On the code before the card 11 of the 12 tests of `tests/r1/test_c34_verify_rule.py`
and `tests/r2/test_c34_stage_labels.py` failed, and both tests of `tests/r2/test_c34_labels_memory.py` failed
(`scripts/pod_pytest.sh C-34 ...`). The first failure:

```text
FAILED tests/r1/test_c34_verify_rule.py::test_c34_a_release_ignores_the_label_files_a_stage_gained_after_it
E   AssertionError: rC34.a
E   assert ['not in mani...rt-0.parquet'] == []
E     Left contains one more item: 'not in manifest: l0seq/contig_label/stage=S1/part-0.parquet'
```

The others failed on `invalid choice: 'labels'` and on `no attribute 'occurrence_files'` or `'write_labels'`. After
the card all pass. What they show:

- **Oracle.** The step's label files have the sha256 of a re-ingest's with ingest's label writer:
  - on card W-36's stages (S1 of `part`; S0 and S1 of `none`; S1 of `part` equals S1 of `full`);
  - on a constructed stage S2. It holds T4, a contig S0 owns (DEC-26); a planted lambda in three records, one object
    listed under two sources; and parts cut one row each.
- **Determinism.** The label and ledger files are byte-identical at 1 and 4 DuckDB threads, with the occurrence files
  listed in reverse, and in a child process.
- **Refusal.** The step refuses (exit 2, nothing written) a stage labelled by ingest, a stage already labelled by the
  step, and a stage with no `contig_occurrence`.
- **Run record.** Under `run-step` the step publishes the label file and the ledger with a run record. A second call
  resumes; the ledger balances (L-01).
- **Readers of a published release** (`rW36.a` on `part` and `none`). `find` of a contig of each stage, `seq
  --contig`, `search hits` TSV and FASTA, the catalog's `contig_label` view, `release verify` (exit, stdout,
  stderr) and `releases` (text and JSON) give the same bytes before and after the stages gain label files.
- **A release defined afterwards** over the same stages lists every stage's label file. It verifies with no label
  warning, `releases` shows `contig_label yes` (the published release still shows `partial(not S1)` or `no`), and
  its label view has the rows of the published release's.
- **Memory** (AGENTS §4.2; `tests/r2/test_c34_labels_memory.py`). The stated bound is 1.25 GiB (*derived*: DuckDB's
  1 GB limit plus 256 MiB). Peak RSS of the step in a child process, by `runstep.step_peak`:
  - 155.3 MiB on the 100x synthetic store (its label file equals ingest's by sha256);
  - 503.0 MiB on a constructed stage of 1,000,000 records over 600,000 contigs.

## Rehearsal (scratch root of links, nothing real written)

`/tmp/if2453/c34/rehearsal` holds links to each stage's `contig_occurrence` and ingest ledger, its own `runs/` and
`ledger/stage=<S>/`, and a scratch profile. The step command is the real run's with that root, units `c34-rh-<s>`,
03:03 to 03:05:

```text
pod_run.sh --profile <scratch> --image bio-box_faec1f4.sif --digest sha256:d67300aa... --mem 8G --threads 2 \
  --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/c34/tree/src --log ... c34-rh-s1 run-step --profile <scratch> \
  --root /tmp/if2453/c34/rehearsal --stage S1 --step labels --threads 2 -- python -m biobox.cli stage labels \
  --root /tmp/if2453/c34/rehearsal --stage S1 --cpus 2 --tmp /tmp/if2453/c34/scratch-rh/S1.labels
```

| Stage | Records | Contigs (label rows) | Objects | Label file | Ledger | Step: wall, CPU-min, peak RSS | Unit: CPU-s, MemoryPeak |
|---|---|---|---|---|---|---|---|
| S1 | 34,076 | 27,258 | 1 | 1,224,492 B | 2,196 B | 3.11 s, 0.054, 324.8 MiB | 8.35, 611 MB |
| P1 | 539 | 539 | 250 | 34,774 B | 10,911 B | 2.85 s, 0.049, 324.5 MiB | 7.58, 609 MB |
| P2 | 5,018 | 5,018 | 2,247 | 311,829 B | 84,781 B | 3.27 s, 0.056, 324.9 MiB | 8.60, 609 MB |

- The contigs of the three stages are 32,815, the number card W-36 *derived* from `rN.3`.
- The step's figures are from its run records. The unit's figures are `pod_run`'s usage line, which counts run-step's
  start-up and its 5 s CPU reading; each unit took about 23 s of wall time.

## Oracle on S1 (a re-ingest into scratch)

The re-ingest read S1's inputs, the object its ingest run record names (`--objects
/tmp/if2453/r10/inputs/S1/objects.json`: `phage_first34076.fa`, 2,121,124,779 bytes, read-only). It ran `bio-box stage
ingest` with the card's code into `/tmp/if2453/c34/oracle/root`, in unit `c34-oracle-s1`: 03:18:17 to 03:25:56, 7.15
CPU-minutes (429.14 CPU-s), MemoryPeak 3.41 GB. S1's original ingest took 7.54 (its run record).

`c34_oracle.py` (`/tmp/if2453/c34/oracle/oracle.out`):

- **The card's oracle.** The label file of `stage labels` on the root's S1 and the re-ingest's label file have the
  same sha256, `a297d9831c3ee504f0ea7af7b1a9e9dfc5db623525b0a98bb44dc14e94581aed`: 27,258 rows each, equal.
- **The same inputs.** The re-ingest's `contig_occurrence` equals the root's S1 file in rows: 34,076 rows, every
  column in schema order, rows in file order, compared as Python values. It also has the same sha256,
  `fcbeef96831f7ddd1ce7aed7a89e2b52151f060152f147d3778357a5dc0a2533`, although the root's file was written in the
  older image. The Parquet writers give the same bytes for this table.
- Last line: `oracle: same`.

## The real root

The driver modes ran from 03:26 to 03:48. Before each launch: no `STOP` in the root, and `top -u if2453` showed at
most about 2 CPUs of the account's other work (this card's test run, card E-16). Every unit was capped (pod_run
`--mem 8G`, the driver's scopes 4 to 16G) at nice 19, with 2 threads. Nothing of the root was deleted or rewritten,
no release was defined in it, and nothing under a `release=` directory was written.

**Before** (`before`, 03:26:56 to 03:33:55).

- `release verify` of the six releases, with the card's code:

  | Release | Stages, files, artifacts | Manifest hash |
  |---|---|---|
  | `rN.1a` | S1, 23 files, 4 artifacts | `39020b27de1fa163570c4255365aebdd` |
  | `rN.1b` | S1,P1, 37 files, 7 artifacts | `d3212c9f264e25208158534e61e24011` |
  | `rN.1` | S1,P1,P2, 52 files, 34 artifacts | `84c2c057b09643bb4e9925938703095d` |
  | `rN.2` | 59 files, 34 artifacts | `2b4b0e43cecdb3a5dbccd1a2bb2b0393` |
  | `rN.3` | 60 files, 34 artifacts | `eee598acd16324d1dd73349cfd6d931f` |
  | `rM.2` | S1,P1,P2,P3,P4,P5,P6, 158 files, 154 artifacts | `2239df7b92ec775058fb9c2ddae3609e` |

  Each line read `verified`, exit 0.
- `bio-box releases`.
- `find --accession NC_001604` and `NZ_AP026781.1` on `rM.2` and `rN.3`: 50 and 1,032 lines each.
- The snapshot: 3,403 files and 80,996,171,076 bytes, with sha256 of the 336 files under a `stage=S1|P1|P2`
  directory.

**Run** (`run`, units `c34-s1`, `c34-p1`, `c34-p2`, 03:35:57 to 03:37:24, one at a time; allowance 55 threads at
launch). Each published 2 files:

| Stage | Step: wall, CPU-min, peak RSS | Unit: CPU-s, MemoryPeak |
|---|---|---|
| S1 | 2.65 s, 0.046, 325.4 MiB | 7.97, 609 MB |
| P1 | 2.44 s, 0.042, 324.0 MiB | 7.25, 609 MB |
| P2 | 2.75 s, 0.048, 324.5 MiB | 7.53, 608 MB |

- The run records name image `sha256:d67300aa…` and `code_id` `a76fd7d3f758…`.
- The six files are byte-identical to the rehearsal's, and S1's label file is the oracle's, `a297d983…`.
- Cost of the three steps together: 0.136 CPU-minutes (*derived*: the sum of the run records' `cpu_min`); of the
  units, 22.7 CPU-seconds (*derived*: the sum of their `CPUUsageNSec`).

**After** (`after` and `check`, 03:38:13 to 03:47).

- The snapshot: 3,412 files and 80,997,844,111 bytes. `c34_check.py` passed (`check: ok`):
  - every one of the 3,403 files is there with the same bytes and mtime;
  - the 336 stage files have the same sha256;
  - the 9 new files are exactly the three label files, the three `labels` ledgers and the three run records.
- **Added: 9 files, 1,673,035 bytes**, listed with sha256 in `2026-11-15_label_files_old_stages/root_added.tsv` and
  on the pod in `/tmp/if2453/c34/root_added.tsv`. They are not yet in the durable copy on `/stor`. A tar of them alone
  would be under 64 MB (AGENTS §7.3), as would card C-33's 56 files, so they go with the next addition to the root.
  They can be made again in seconds from the stages' `contig_occurrence`.
- `release verify` of the six releases: each `verified`, exit 0, with the same file counts and manifest hashes as
  before.
  - Their stdout and stderr equal the before logs except one line. For a release over S1, P1 or P2, the `built by`
    line now counts the 3 `labels` run records of code `a76fd7d3f758`. For `rM.2`: `built by: 6d3c4606bd13 (..., 14
    steps), a76fd7d3f758 (..., 3 steps), ec1272b840c6 (..., 449 steps); 536 steps unknown`.
  - `rM.2` keeps W-36's warning that S1, P1 and P2 have no `contig_label` files: its manifest is as it was.
- `find` of the two controls on `rM.2` and `rN.3`: byte-identical to the before outputs, stdout and stderr (`diff -r`
  of the two directories).
- `bio-box releases`: the published releases show what they showed. `rM.2` is `partial(not S1,P1,P2)`; `rN.1a`,
  `rN.1b`, `rN.1`, `rN.2` and `rN.3` are `no`. Only `rM.1` changed, from `partial(not S1,P1,P2)` to `yes`.

```text
rM.1	no	S1,P1,P2,P3,P4,P5,P6	-	-	45815843	pfam-38.2	-	yes	yes	no
rM.2	yes	S1,P1,P2,P3,P4,P5,P6	39224	54265	45815843	pfam-38.2	-	yes	partial(not S1,P1,P2)	yes
rN.3	yes	S1,P1,P2	29755	32815	10620684	pfam-38.2	-	yes	no	yes
```

**A release defined afterwards, on a scratch view** (`seven`). The view `/tmp/if2453/c34/view7` links to the root
with its own `release_def`, and the release `rC34.v` was defined there over the seven stages. No `release define` ran
in the root. `bio-box releases --root /tmp/if2453/c34/view7`:

- before the run (03:35:27): `rC34.v	no	S1,P1,P2,P3,P4,P5,P6	-	-	-	pfam-38.2	-	no	partial(not S1,P1,P2)	no`;
- after it (03:48): `rC34.v	no	S1,P1,P2,P3,P4,P5,P6	-	-	-	pfam-38.2	-	no	yes	no`.

The real release over the seven stages is card C-33 part 2.

## Seen

- Six releases, not four, report the new files on the code before the card: `rN.1a` and `rN.1b` list S1 (and P1) too.
- After the files exist, `rM.1` (defined, no manifest) reads its stages' directories (`ReleaseFiles` from its
  `release_def`, card B-28) and shows `contig_label yes` in `bio-box releases`, where it showed
  `partial(not S1,P1,P2)`. The published releases do not change: their manifests do not list the files.
- `release verify`'s `built by` line of a release over S1, P1 or P2 now counts the `labels` run records among its
  stages' steps. The line is information only and comes from run records, not from the manifest.
- The step's own cost is about 3 s and 0.05 CPU-minutes per stage, below the 300 CPU-minute budget by a factor of
  6,000. No estimate is made for it in `build-stage`; it is not one of its steps.
