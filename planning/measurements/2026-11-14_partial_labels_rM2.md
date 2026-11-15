# A release whose stages are only partly labelled: lookups by accession on `rM.2` (card W-36)

Card W-36. The release check of `rM.2` stopped at its first step (*measured*, wilkcomp02, 2026-11-14 22:13 to 22:15,
unit `ex-rM2`, `/tmp/if2453/ex-rM2/out/summary.tsv`: 37 invariants missing). Both worked examples find their positive
control by accession, and `bio-box find` refused both:

```text
refused: no contig of release rM.2 has the accession NC_001604 (the accession column of `bio-box search hits --format tsv`)
refused: no contig of release rM.2 has the accession NZ_AP026781.1 (the accession column of `bio-box search hits --format tsv`)
```

(`/tmp/if2453/ex-rM2/out/rnap/logs/find-t7.log`, `.../crispr_cas_neighbourhoods/logs/find-tniq.log`.)

## Cause

`ls /tmp/if2453/r10/root/l0seq/contig_label` shows `stage=P3` to `stage=P6` only, and `contig_occurrence` has all
seven stages (*measured*, wilkcomp02, 2026-11-14). S1, P1 and P2 were ingested before card B-27, which added the
per-stage `contig_label` table at ingest. Every reader of a contig's accession goes through one of two places:

- the export store's label relation (`biobox.search.export.open_store`, relation `lab`). `find`, `seq --contig`,
  `search hmm --positive`, `search hits --format tsv|faa`, `nbr flanks` and the locus tables read it through
  `_annotate_contigs`. Card W-18's fallback to `contig_occurrence` was taken only when the release listed no label
  file at all (`rN.3`). With label files for some stages, the other stages' contigs had no row, and so no accession;
- the catalog's `contig_label` view (`catalog.label_view_sql`, stored by `generate_catalog`), with the same rule.

`release verify` warned only when a release had no label file at all, and `bio-box releases` printed `contig_label yes`.

## Fix (readers only)

No file of any release is written or rewritten, no manifest or schema changes, and there is no new table.

- `export.open_store`: when the release lists label files for some of its stages and not for others
  (`export._unlabelled`), the label relation is the label files of the labelled stages `UNION ALL` the
  `contig_occurrence` rows of the others, one row per record with `n_records` 1. This is W-18's rule, applied per
  stage. `_annotate_contigs` already merges a contig's rows over stages (first record by
  (`l0_object_sha256`, `record_ordinal`), records summed), so the answer is the one a fully labelled release gives.
  One `note:` line per command names the stages.
- `catalog.label_view_sql(files, n_stages, fallback)`: the same union in the stored view (`catalog.unlabelled_stages`,
  `_label_fallback`). `open_catalog` shadows the stored view of a catalog built before the card with a temporary view
  of the same text (`_shadow_partial_labels`). `rM.2`'s catalog is such a catalog and is not rewritten.
- `release verify`: `warning: release R: stages S1,P1,P2 have no contig_label files (...)`. The exit code is unchanged.
- `bio-box releases`: `contig_label` shows `partial(not S1,P1,P2)`. With `--json` the record gains `unlabelled_stages`;
  `lookups.contig_label` stays `true` when any stage has label files.

On a release whose stages are all labelled, or none is, every path is the one before the card: `_unlabelled` returns
None, `fallback` is empty, and nothing is shadowed.

## Tests (`tests/r2/test_w36_partial_labels.py`)

Three fixture roots are built from the same inputs: S0 (micro-S0 lambda and T4) and S1 (a lambda genome with one base
changed, given by two source objects so its contig has two records), gene-called and taped, with release `rW36.a`
over both, a manifest and a catalog.

- `full`: both stages labelled.
- `part`: S1 ingested with the label writer left out, as ingest worked before card B-27. Nothing was deleted.
- `none`: neither stage labelled.

S0's label file is written one row per row group. The oracles are loops over each stage's own `contig_occurrence`
files and the FASTA written by the test.

On the code before the fix, 7 of the 10 tests failed (*measured*, wilkcomp02, 2026-11-14, `scripts/pod_pytest.sh W-36
tests/r2/test_w36_partial_labels.py`). The first failure:

```text
FAILED tests/r2/test_w36_partial_labels.py::test_w36_find_by_accession_finds_the_genes_of_every_stage[S1]
E   AssertionError: ('part', 'refused: no contig of release rW36.a has the accession lambda_w36_again (the accession
    column of `bio-box search hits --format tsv`)')
    assert 2 == 0
```

The others that failed: `seq --contig ACC:..` (refused), the hit TSV (empty `accession`, `n_records`), the catalog view
(no S1 row), a catalog built by the old code, `release verify` (no warning) and `releases` (no `partial`). After the
fix all 10 pass. In each test the bytes of `find`, `seq`, `hits.tsv`, `hits.faa` and the view's rows are equal on the
`part`, `full` and `none` roots.

**Pruning** (AGENTS §4.4). In `part`, one S0 contig looked up by contig_uid from a joined table (the exports) and as a
literal (the catalog view) reads 1 of the 2 row groups of S0's `contig_label` and nothing of S0's `contig_occurrence`.
The unlabelled S1 is read from its `contig_occurrence` whole, which is W-18's cost for those stages.

On `rM.2` (*measured*, wilkcomp02, 2026-11-14, `/tmp/if2453/w36/bin/prune.py`, in the image with the tree's `src`), one
lookup by contig_uid through the export store reads:

- `contig_occurrence` of S1, P1 and P2: 3 of 3 row groups. Each stage's file is one row group: 32,815 contigs in all
  (*derived*: `rN.3`'s contigs, the same three stages).
- `contig_label` of P3 to P6: 2 of 4 row groups for a P3 contig; 4 of 4 for an S1 or P2 contig, because each stage's
  label file is one row group whose uid range spans nearly all uids.

Each lookup took 0.02 to 0.08 s.

## On the pod (read-only on `/tmp/if2453/r10/root`)

Run with the runtime image `bio-box_faec1f4.sif` (sha256 `d67300aa…`) and the tree's `src` (commit `a123da1`), by
`/tmp/if2453/w36/bin/lookups.sh`, on wilkcomp02, 2026-11-14 (*measured*):

- `bio-box find --release rM.2 --accession NC_001604`: exit 0, 49 genes of contig `719747a8fef97b2380955d92845861db`
  (T7, stage S1), the first `NC_001604:924-1278 +1`; 5.1 s of wall time.
- `bio-box find --release rM.2 --accession NZ_AP026781.1`: exit 0, 1,031 genes of contig
  `9acd61e52f51368dab7340cb1d16544b` (stage P2); 4.5 s.

  Both print `note: release rM.2: stages S1,P1,P2 have no contig_label files (...)`.
- `bio-box releases`: `rM.2 yes S1,P1,P2,P3,P4,P5,P6 39224 54265 45815843 pfam-38.2 - yes partial(not S1,P1,P2) yes
  26e10f05401faeefb78769e83f7f6cb4`. `rM.1`, clustered only, shows `partial(not S1,P1,P2)` too.
- `release verify --release rM.2`: `release rM.2: verified (stages S1,P1,P2,P3,P4,P5,P6, 158 files, 154 artifacts,
  manifest hash 2239df7b92ec775058fb9c2ddae3609e)`, exit 0, with `warning: release rM.2: stages S1,P1,P2 have no
  contig_label files (ingested before card B-27); readers take those stages' labels from their contig_occurrence,
  which cannot prune a lookup by contig there`.
- `release verify --release rN.3`: `verified (... manifest hash eee598acd16324d1dd73349cfd6d931f)`, exit 0, with
  W-17's warning as before.

**Byte identity on `rN.3`** (a release with no label file). Both worked examples were run on `rN.3` twice at once,
with the runtime image above (*measured*, wilkcomp02, 2026-11-14 22:38 to 23:10). Each run is
`scripts/release_examples.sh rN.3` in a unit of its own, 2 threads per launch:

- `w36-ex-rN3-new`: this card's tree, `/tmp/if2453/w36/tree`.
- `w36-ex-rN3-base`: the tree before it, `66ef9c7`, unpacked from `git archive` into `/tmp/if2453/w36/base`.

`/tmp/if2453/w36/bin/compare.py` hashed every file of each example's output, leaving out the logs, scratch and
search-run directories.

| Example | Files | Identical | Differing |
|---|---|---|---|
| CRISPR-Cas | 3,567 and 3,567 | 3,559 | 8 |
| RNA polymerase | 5,628 and 5,628 | 5,623 | 5 |

- The 13 differing files are the 7 `RUN` files (the run id), 4 `recipe.yaml` and 2 `hmm/SHA256SUMS` (the path of
  `OUT`), and 2 `oracle.tbl` (hmmsearch's header: paths and date), all read by `diff`.
- Every table, figure, locus file and SQL output is identical.
- `summary.tsv` is identical: 92 checks, every invariant holds in both runs.
- CPU per run, CRISPR-Cas and RNA polymerase (*derived* by the script from the units' usage lines): new code 1,355.1
  and 1,114.6 CPU-seconds; old code 1,379.6 and 1,109.2.

**The release check of `rM.2`** with this code (unit `w36-ex-rM2`, 22:38 to 00:06) ended `release_examples.sh: rM.2:
2 examples, every invariant holds`: 37 invariants hold and 55 recorded values are new
(`2026-11-14_release_examples_rM2.md`).

## Proposal (step 5 of the card, not built): label files for old stages without a re-ingest

1. **A stage step `bio-box stage labels --stage S`**, one `run-step` per stage.
   - It reads the stage's `contig_occurrence` files listed in its ingest run record (an explicit list, AGENTS §3.9).
   - It sorts the rows by (`contig_uid`, `l0_object_sha256`, `record_ordinal`) in DuckDB, spilling to scratch.
   - It writes `l0seq/contig_label/stage=S/part-<k>.parquet` through the function ingest uses (`store._label_batches`
     and `io.write_parts`): new files only, under a temporary name and renamed, with a ledger row.
   - Its oracle: re-ingest a stage into scratch from the same objects, and its label files must have the same sha256
     as the step's.
   - Cost: one sort of 32,815 contigs' records for S1, P1 and P2 together, seconds (*estimate*).
2. **A rule in `verify_manifest`.** Today a new `contig_label` file in S1 would be "not in manifest" for every release
   that lists S1 (`rN.1`, `rN.2`, `rN.3`, `rM.2`). `contig_label` has no partition but its stage, so card C-19's
   `_owner_key` makes it every listing release's.
   - The rule: an unlisted `contig_label` file of a stage for which the manifest lists none is not the release's
     concern. That release was built before the stage had labels, and its readers take that stage's labels from
     `contig_occurrence` (this card).
   - Without that rule the step must not run on a stage that a verified release lists.
3. **Then the release after `rM.2`** (card C-33) is fully labelled, and its lookups prune on every stage. Its readers
   give the same bytes as on `rM.2`: the tests here show a partly labelled release and a fully labelled one agree.
