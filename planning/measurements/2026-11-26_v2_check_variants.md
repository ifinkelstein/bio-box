# V2's check and the non-default variants of contigs another stage stores (card R-17, 2026-11-26)

`V2`'s `check` stopped wave-4 at 16:32:54 CST with "contigs both stored by the stage and held by another: 43"
(`/tmp/if2453/stages/v2/check/1/check.out`; `planning/tasks/STATUS.md` W10). This note records why that was the check's
fault and not the stage's, what the check asserts now, and the changed check's read-only run on `V2`.

## 1. What a later stage's non-default variant is (card step 1)

The plan documents and the code define it; no decision row or test says otherwise.

- `planning/data_layer_and_tests.md` §2A.3: "`contig_index` maps `contig_uid → raw_digest`. Where several raw variants
  share a `contig_uid`, the default variant is the one with the lowest (source precedence, `l0_object_sha256`,
  `record_ordinal`), so it does not depend on input order (MR-01), and the other variants are reachable by
  `raw_digest`. `contig_index` is partitioned by `stage`; a later stage adds only new (contig_uid, raw_digest) pairs
  and never changes an earlier default (implemented in `src/biobox/l0seq/store.py`, tested by
  `tests/r2/test_micro_s0_l0seq.py`)." The paragraph before it names the case: "a soft-masked copy and an uppercase
  copy, MR-14".
- `src/biobox/registry/schemas.py`, `contig_index`: primary key and sort key `(contig_uid, raw_digest)`, partitioned by
  `stage`; doc "Default variant: lowest (source precedence, l0_object_sha256, record_ordinal)".
- `src/biobox/ids.py`, `contig_uid`: "BLAKE3-128 of the canonical uppercase sequence." So two raw sequences with one
  `contig_uid` have the same bytes in upper case: they differ in letter case only, at the same length (*derived*).
  A rotated circular record is another `contig_uid`, not a variant; card R-17's "as a rotated circular record gives
  it" does not hold, and the new test uses a soft-masked copy.
- `src/biobox/l0seq/store.py`, `_write_tables`: "contig_index: new (contig_uid, raw_digest) pairs; the default is the
  variant with the lowest (precedence, object sha256, record ordinal) of the stage, unless an earlier stage gave
  one"; its SQL writes a pair no earlier stage has (`NOT EXISTS ... old_idx`) with `is_default_variant` false when an
  earlier stage holds the contig's default (`old_default`).
- `store.ContigOwners` (the owner of DEC-26): "the stage whose `contig_index` holds the contig's default variant.
  Ingest writes a sequence once, in the first stage that sees it, and a later stage's `contig_index` holds no default
  variant of a contig an earlier stage stored, so that stage is unique". Gene calls, IGRs, tapes and CRISPR skips all
  name that owner, so a later stage's non-default row owns nothing.
- DEC-26 ("A contig shared by genomes in different stages belongs to the stage that stored it first") is about
  ownership. It does not stop a later stage from storing new raw bytes of the contig.

Tests:

- `tests/r2/test_micro_s0_l0seq.py::test_mr14_soft_mask_keeps_contig_uid`: one stage, lambda and a soft-masked copy
  from a lower-precedence source. Same `contig_uid`, other `raw_digest`, the masked bytes readable by `raw_digest`,
  and the default is the higher-precedence variant.
- `tests/r2/test_micro_s0_catalog.py::test_mr14_macro_reads_only_the_default_variant`: the SQL macros read the
  default, and `seq_raw` reads the variant.
- `tests/r2/test_micro_s0_contracts.py::test_r09_a_later_stage_keeps_earlier_files`: stage `S1` ingests a soft-masked
  lambda that `S0` stores ("the lambda variant is a new S1 shard"), and the append-only check passes. This is the
  case `V2` met. No test before card R-17 asserted the flag of that row. The new fixture test does, from ingest's own
  output: `tests/r1/test_r17_stage_check_variants.py::test_the_fixture_is_a_later_stage_holding_a_new_raw_variant_of_a_contig_the_earlier_stage_stores`.

## 2. The ownership block of `stage_check.py`, before and after

Before the change, `ci` held every contig with any row in the stage's index. Now `ci_all` is that set and `ci` holds
only the contigs whose default variant the stage has. `held` is unchanged: contigs of `genome_input` whose default
variant another stage has.

| Assertion (its failure line) | Before | After |
|---|---|---|
| contigs held by two other stages (DEC-26) | `held` with 2+ stages | unchanged |
| contigs both stored by the stage and held by another | any row of the stage × `held` | a default-variant row of the stage × `held`; same message |
| contigs of the stage's index without exactly one default variant among all stages (DEC-26) | none | new: for each contig of the stage's index, its default rows here plus in every other stage = 1 |
| rows of the stage's index whose (contig_uid, raw_digest) another stage's index holds | none | new: every row, default or not (`(contig_uid, raw_digest)` is the table's key across stages; ingest writes only new pairs) |
| `<n>` non-default rows in the stage's index, but the index files of `[...]` hold no raw_digest | none | new: a non-default row whose digest cannot be compared fails (only hand-built roots lack the column) |
| genome_input contigs neither stored nor held | `gi` − any row − `held` | `gi` − default rows − `held` (stricter) |
| stored contigs not in genome_input | every contig of the index | unchanged (`ci_all`) |
| ledger skips for another stage that differ from the stage holding the contig | the whole reason against `held` | unchanged |
| stored contigs under the least length | every contig of the index | unchanged (`ci_all`) |
| coverage (genes or a ledger row; no genes on a contig not stored) | every contig of the index | default-variant contigs (`ci`) |

Summary keys: `stored` is now the default-variant count, the count the estimates use (`stage.stage_contigs`,
`bases_stored`). New keys beside it: `index_contigs`, `index_rows`, `non_default_rows_of_stored_contigs`,
`variants_of_contigs_held_by_other_stages` (`contigs`, `rows`, `bases`, `by_stage`), and one count per new
assertion. No key was removed or renamed.

Readers of the report (`rg` over `scripts/`, `src/`, `tests/`, `RUNBOOK.md`, 2026-11-26). `stage_driver.sh check`
reads the exit status, the last line of `check.out` and `root_added.tsv`. `stage_queue.sh`, `stage_wave.sh`,
`stage_est.py`, `release_train.sh` and `train_check.py` read nothing of it. `tests/r1/test_stage_driver.py` reads
`failed`, `contigs.ledger_skipped_other`, `contigs.ledger_skipped_within` and `pfam.copy_check`.
`tests/r1/test_a21_stage_check_copy_sources.py` reads `failed` and `pfam`. No code reads `contigs.stored`.

## 3. The changed check on `V2`, read-only

*Measured* (`wilkcomp02`, 2026-11-26, 17:40:00 to 17:41:43 CST, about 102 s from the scratch directory's creation to
`summary.json`). The command `stage_driver.sh check` runs (`bb py stage_check.py ...`, the same image, binds and
scope, 2 threads) ran on the card's branch as mirrored by `pod_pytest.sh` to `/tmp/if2453/ci/R-17`, with every
output and the profile's `tmp` under `/tmp/if2453/r17/`:

    systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
      singularity exec --cleanenv --env PYTHONPATH=/tmp/if2453/ci/R-17/src --env OMP_NUM_THREADS=2 \
      --env PYTHONUNBUFFERED=1 --bind /tmp/if2453 --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro \
      --bind /sys/fs/cgroup:/sys/fs/cgroup:ro /stor/home/if2453/bio-box/images/bio-box_faec1f4.sif \
      python /tmp/if2453/ci/R-17/scripts/stress/build_stage_pod/stage_check.py --conf /tmp/if2453/r17/V2.conf \
      --stage V2 --root /tmp/if2453/r10/root --before /tmp/if2453/r17/snap/before --after /tmp/if2453/r17/snap/after \
      --objects /tmp/if2453/s2/inputs/V2/objects.json --out /tmp/if2453/r17/out --tmp /tmp/if2453/r17/tmp/check \
      --profile /tmp/if2453/r17/profile.yaml --before-est /tmp/if2453/r17/snap/estimate/before.json

Inputs: `V2.conf` copied from the tree the wave ran (`/tmp/if2453/ci/rel1126a`), sha256 `63d4d061...ee1b`, equal to
the conf hash in `after/DONE`. `listing.tsv` and `earlier.sha256` of `before` and `after`, and
`estimate/before.json`, were copied with `cp -p`, and their sha256 equal the originals. The profile is `V2`'s
`profile.yaml` with `tmp: /tmp/if2453/r17/tmp`. Afterwards, `find /tmp/if2453/stages/v2 -newer <the scratch
profile>` and `find /tmp/if2453/r10/root -newer <the same> -print -quit` printed nothing.

Result: `RESULT: all checks pass`. Counts (`/tmp/if2453/r17/out/summary.json`):

| | count |
|---|---|
| `genome_input` rows, contigs, genomes | 261,578; 261,578; 230,050 |
| stored (default-variant contigs) | 243,466 (card's expectation 243,466) |
| bases stored | 8,199,089,073 |
| index contigs, index rows | 243,509; 243,915 |
| non-default rows of contigs `V2` stores | 406 |
| variants of contigs held by other stages | 43 contigs, 43 rows, 320,799 bp, all of `S1` (card: 43 of `S1`) |
| held by other stages | 18,112 (card: 18,112): `V1` 16,803, `H1` 1,136, `S1` 70, `PL1` 68, `P4` 10, `G1` 8, `SPX1` 7, `P3` 5, `P2` 3, `P5` 2 |
| ledger skips for another stage | 18,112 |
| held by two other stages; both stored and held; not exactly one default; a repeated pair; neither stored nor held; stored not in `genome_input`; ledger skips that differ | 0 each |
| new files, all under `stage=V2` | 194 (6,955,837,434 bytes) |

Check: 243,466 + 43 = 243,509 index contigs, and 243,466 + 406 + 43 = 243,915 rows (*derived*). Both agree with
the diagnosis in STATUS W10.

## 4. The checks written before `stage_check.py`

`c47_check.py` has the same rule: `owned` is every contig of the stage's index (line 284), and `owned & stored_by`
fails (line 298). `c48c_v1_check.py` (lines 253, 273), `c26_check.py` (250, 266) and `c24_check.py` (227, 245) have it
too. They are run only by their own drivers (`c47_driver.sh check` and `rehearse`, `c47_queue.sh`; `c48c_v1_driver.sh`,
`c26_driver.sh`, `c24_driver.sh`), which RUNBOOK ("Building a stage") keeps only for the stages they built. No current
script, queue or wave calls them. They were not changed (not in the card's Files).
