# S2C phase 2: skani pinned, ANI selection code, the genome fetch (card C-27 steps 1 to 3; 2026-11-12/13)

Steps (1) to (3) of card C-27 (DEC-31). Step (4), running skani, waits for an image built with the new pin; it is
prepared (driver, dry-run plan) and was not run. Every number is labelled *measured* (host, date, command),
*derived* (formula and inputs) or *estimate*. Times are CST. Nothing was written under `/stor`, `/tmp/if2453/r10`
or `/tmp/if2453/r50`.

## Result in one paragraph

skani **0.4.0** (bioconda, linux-64 build `hec9b1f2_0`) is pinned. `c25_select.py --ani` picks 20 per species by
greedy farthest-first from the best-quality genome. It matches a separately written brute-force oracle and does
not depend on input order, and without `--ani` its output is byte-identical to before. The skani driver
`c27_skani.py` is written and tested against a stand-in binary, and its dry run on the real candidates plans
499,624 farthest-first pairs in 5 shards (about 1,048 CPU-minutes at the upper cost *estimate*). **The fetch is
done** (2026-11-12 23:17 CST; a chain of detached units, one per batch of about 9.9 Gbp): 28,457 of the 28,483
genomes asked were served and verified (118,675,230,482 bp, 34,788,621,074 bytes kept as `.fna.gz`), 26 were not
served, and the 7 phase-1 retries were all refused again, so they stay out per DEC-31. No refusal, no md5
mismatch. The orchestrating session filled in the per-batch table below from each batch's `fetch_result.json`.

## Step 1: the pin

`container/environment.yml` gains one line: `skani=0.4.0`, with a comment naming DEC-31 and card C-27.

- Source: bioconda's package metadata, `https://api.anaconda.org/package/bioconda/skani` (one `curl -sS` GET from
  the laptop, 2026-11-12 about 21:05 CST, curl's default user agent, nothing about us; 66,745 bytes; *measured*).
- Versions listed: 0.1.5, 0.2.0, 0.2.1, 0.2.2, 0.3.0, 0.3.1, 0.3.2, **0.4.0** (`latest_version` 0.4.0).
- 0.4.0 for linux-64: `linux-64/skani-0.4.0-hec9b1f2_0.conda`, build `hec9b1f2_0`, build number 0, 759,597 bytes,
  md5 `d229015778df0c1eec388ee9c974a6f8`, sha256 `7c2041a84c20ca903a363f032099dbc9c4443490a84f3b52c1915d15760c326c`,
  uploaded 2026-10-11 06:08 UTC (a month before the pin, so older than the two-week rule for pins); depends
  `__glibc >=2.17,<3.0.a0`, `libgcc >=14`, `libstdcxx >=14`; licence MIT. Builds exist also for linux-aarch64,
  osx-64 and osx-arm64.
- Why this one: the card asks for the current bioconda release. GTDB R226 drew its species with skani 0.2.1 and
  GTDB's web service runs 0.3.0 (`2026-11-12_s2c_ani.md`); nothing here needs their exact values, only a
  within-species ordering, so the newest release is pinned.
- skani was not installed or run anywhere. It runs only from the image rebuilt when `main` next takes
  `plan/task-cards` (AGENTS §7.1).

## Step 2: the selection code

### `c25_select.py --ani TABLE`

- Species with more than 20 complete genomes: 20 picked by greedy farthest-first from rank 1 (the best-quality
  genome; so every phase-1 genome stays selected, and the script asserts it: every rank-1 genome is picked, and
  every phase-1 genome the quality order selects is in the ANI selection).
- Distance of a pair: 100 − ANI. After each pick every unpicked genome keeps its distance to the nearest pick; the
  next pick is the genome farthest from all picks. Ties by quality order (CheckM2 score, then GTDB accession).
- **A pair skani does not report counts as the largest distance (100)**: a row with an empty ANI (skani was asked
  and printed nothing: below its ANI screen), a pair absent from the table (with `--require-pairs`, an absent
  consulted pair is an error instead), or a row whose larger alignment fraction is under `--min-af`.
- **Low alignment fraction.** `--min-af` defaults to 15 and is applied to the larger of the two fractions: an ANI
  computed over less than 15% of both genomes is not a genome-wide value, so it is treated as not reported. 15 is
  skani's own default for printing a value at all (its `--min-af`, as its documentation describes it: printed when
  one genome's aligned fraction exceeds it; to be checked against `skani dist --help` in the image), so with
  skani's defaults the rule only matters for a table made with a lower setting. Above 15 the fraction is not used:
  the plan's rule is ANI diversity (`STAGES.md`), and in complete genomes of one species a lower fraction mostly
  measures accessory content (plasmids, islands), not the core divergence the rule is meant to spread. The
  fractions stay in the table so step (4) can report them for the picked pairs.
- A pair given twice (for example in both directions) keeps its larger ANI, so row order cannot matter.
- Outputs, with `--ani`: `s2c_selected.tsv` and `remainder.tsv` gain `ani_pick` and `ani_to_nearest_earlier_pick`;
  `ani_species.tsv` gives per species how many picks differ from the quality order and the minimum pairwise ANI
  among the 20 under each rule (with how many of its 190 pairs are reported, unreported or absent);
  `selection.json` gains an `ani` block (consulted pairs, genomes that differ, phase-1 subset, quantiles of the
  per-species minimum ANI under each rule, the ten largest species).
- Without `--ani` the outputs are byte-identical to the script before this card (test below).

### `c25_select.py --candidates`

Writes `ani_candidates.tsv` (every complete genome of the species over the cap) and `ani_fetch.tsv` (those not in
phase 1, plus the `--phase1-not-served` genomes with `fetch_reason = phase1_not_served_retry`).

### `c27_skani.py` (the driver for step 4; written, tested against a stand-in binary, not run for real)

- **Pairs: the farthest-first rounds, not all pairs.** 19 sequential rounds per species: `skani dist` with the
  newest pick as query against the sketches of the genomes not yet picked; the next pick is chosen by the same
  function `c25_select.py --ani` uses (`farthest_first`, reading rows through the same `ani_value`). That is
  Σ (19 n − 190) = **499,624 pairs** (*derived*, `ani_needs.json` per-species counts; the dry-run below recounts
  it), against 15,835,149 for all pairs. Plus one `skani triangle -E` of each species' 20 best-quality genomes
  (at most 190 pairs per species, 46,170 in all, *derived*), for the comparison of the minimum pairwise ANI under
  both rules; its pairs are added only where the rounds did not already hold them, so the table's values for
  the pairs the pick consulted are exactly the ones the driver used.
- One species at a time: its sketches on disk under `WORK/species/<id>/`, its pairs in memory (at most 75,943 for
  *E. coli*, *derived*). Each finished species writes `pairs.tsv`, `species.json` (wall times of every skani call)
  and `DONE`, so a shard resumes.
- Shards: species packed by first-fit decreasing on an estimated cost, (genomes × sketch CPU-s + pairs × pair ms)
  / 60, under a budget of 240 CPU-minutes per shard process (the rule is 300); `plan` refuses a species estimated
  at 300 or more. Units `c27-skani-<k>` (systemd-run --user, nice 19, MemoryMax 32G, CPUQuota threads × 100%,
  inside the image via `singularity exec`); `launch --parallel 2` starts at most two at a time (2 × 8 threads).
  Markers: `WORK/shards/<k>/DONE` or `FAILED`; then `collect` writes the one table for
  `c25_select.py --ani TABLE --require-pairs`.
- skani's command lines and output columns are taken from its README, not from the binary (not available here).
  Before the first real run, step (4) must compare them with `skani sketch|dist|triangle --help` in the new image.

### Dry run on the real candidates (*measured*, wilkcomp02, 2026-11-12 21:34 CST)

`python3 -I c27_skani.py plan sel/ani_candidates.tsv sel/paths_now.tsv /tmp/if2453/c27/skani --image
IMAGE_WITH_SKANI --dry-run` (output `/tmp/if2453/c27/logs/skani_dryrun.txt`; the paths came from the phase-1
locks and B01's lock, so 26,087 candidates had no file yet):

| Count | Value |
|---|---|
| species, genomes, bp | 243; 28,726; 119,700,310,300 (equal to `ani_needs.json`) |
| farthest-first pairs | **499,624** (equal to C-25 1b's derived count) |
| quality-top pairs (triangle) | 46,170 |
| all pairs (not computed) | 15,835,149 |
| estimated CPU-minutes (2 CPU-s per genome, 10 ms per pair: the upper end of C-25's *estimate*) | 1,048.5 (*estimate*) |
| shards under the 240 budget | 5: 239.5, 240.0, 240.0, 240.0, 89.1 CPU-min; 3, 8, 30, 109, 93 species |

Before the real run, step (4) must: build the plan for real with `--image` set to the new image, after every
batch is DONE (plan refuses missing files); measure one species first (`run-shard` on the shard that holds *S.
aureus*, or a one-species plan) against the estimate; then `launch --parallel 2`.

### Tests

`tests/r1/test_c25_select.py` (covers D-01, D-04; ring R1), 13 new tests:

- `test_ani_farthest_first_matches_brute_force_oracle`: an ANI world built in the test (5 species: 14, 9 and 6
  genomes over a cap of 5, one at 5, one of 2; ANI drawn from six two-decimal values so ties occur; 8% of pairs
  unreported, 4% with both fractions under 15) gives the picks of an oracle written separately (quality rank by
  the pairwise `_beats` count; exact `Fraction` distances; at each step every candidate's distance to the picks
  recomputed from scratch). Also: the rule differs from quality order in this world; phase 1 a subset; consulted
  pairs equal Σ ((cap − 1) n − cap (cap − 1) / 2); species at or under the cap unchanged; the remainder; the
  per-species minimum pairwise ANI under both rules against the world.
- `test_ani_min_af_changes_the_rule_as_stated` (`--min-af 0`); `test_ani_input_order_does_not_change_output`
  (metadata rows shuffled, table rows shuffled and a/b swapped: byte-identical outputs);
  `test_ani_absent_pairs_count_as_largest_or_refuse`; `test_ani_table_with_a_stranger_is_refused` (a pair across
  two species, an accession not in the metadata).
- `test_quality_only_output_unchanged_without_ani`: sha256 of the four outputs on the step-1 fixture equal those of
  the script before this card (commit `00e313e`'s file, run on the pod before the change). A no-change check, not
  an oracle; correctness of those outputs is the brute-force test of step 1.
- `test_candidates_and_fetch_list`.
- `test_skani_normalise_hand_made_file`: a file in skani's documented `dist` format (names with spaces, a reversed
  duplicate, a `.sketch` path, an extra column) normalises to the expected table; `--expect-*` fills unreported
  pairs; a file without the fraction columns is refused. `test_skani_normalised_world_selects_like_oracle`:
  skani-format text of a world → normaliser → `--ani` gives the oracle's picks.
- `test_skani_driver_dry_run_plan`: the dry run writes nothing; species, genomes, farthest-first pairs (by
  simulating the rounds), quality-top pairs and all pairs; shards under the budget; the cost formula; the unit
  commands; `launch --dry-run` prints two.
- `test_skani_driver_run_with_stand_in_binary`: `plan`, `run-shard`, `collect` against a stand-in `skani` (a test
  script printing a known world in skani's format): `c25_select.py --ani --require-pairs` on the collected table
  gives the oracle's picks; one `sketch`, 4 `dist` and one `triangle` per species; a second run asks nothing.
- `test_fetch_keeps_members_compressed_and_verified` (a hand-made Datasets-style zip: deflated and stored members,
  a bad md5), `test_fetch_batch_retries_lists_failures_and_resumes` (a stand-in fetcher: batch, single retry, md5
  refetch, failure list, resume, STOP), `test_fetch_split_into_batches`.

They failed first: with the scripts unchanged, `sh scripts/pod_pytest.sh C-27 tests/r1/test_c25_select.py` gave
`13 failed, 5 passed` (`c25_select.py: error: unrecognized arguments: --ani ...`, `... --candidates`, `c27_skani.py`
exit 2, `No such file or directory: .../c27_fetch.py`). After: `sh scripts/pod_pytest.sh C-27
tests/r1/test_c25_select.py tests/r0` gave `213 passed in 10.51s` (wilkcomp02, in the image). `sh -n` and `dash -n`
pass on `c27_fetch.sh`.

## Step 3: the fetch

### The list

`c25_select.py --candidates --phase1-not-served <the 7>` (unit `c27-select`, 1 CPU, 31.4 s, peak RSS
253,568 kB, *measured*). Its inputs were the GTDB metadata C-25 unpacked (sha256 values equal to `sources.yaml`)
and phase 1's `all.tsv` (`e68d2576…`). `s2c_selected.tsv`, `remainder.tsv` and `by_phylum.tsv` are byte-identical
to C-25's (`9a67b230…`, `7ae727a8…`, `3b3e831c…`).

- `ani_candidates.tsv`: 28,726 genomes, 119,700,310,300 bp, 243 species; 243 held by phase 1 (sha256 `d49ee8b9…`).
- `ani_fetch.tsv`: 28,490 rows, 118,833,937,433 bp, sha256
  `74f5429a45be0244729cda063a922bf3e7a1e611c99759ca5ad71f2538290843`. That is 28,483 candidates (equal to C-25 1b)
  plus the 7 retries.
- Split (`c27_fetch.py split --batch-gbp 10`): B00 holds the 7 retries (30,138,452 bp). B01 to B12 hold 2,240
  to 2,639 genomes and 9,896,810,610 to 9,905,019,239 bp each (`/tmp/if2453/c27/batches/split.json`).

### How

`c27_fetch.py` uses phase 1's source, client, pace and verification (it imports `r50_fetch.py`'s `Fetcher`, URL,
STOP path and limits): the `ncbi_datasets_genome` URL from `sources.yaml`, batches of 20, at most 3 requests per
second, the same generic user agent, md5 against `md5sum.txt`, and one retry alone. The FASTA is kept as
`.fna.gz`: each zip member's deflate stream is copied into a gzip container, with no recompression. The file is
then read back through gzip, which checks the CRC, and its md5 is checked. `lock.json` keeps phase 1's fields
plus `gz_bytes` and `gz_sha256`. Offline check: `c27_fetch.py selftest` on the ladder's two-genome zip matched
phase 1's md5 values (`1fd428ef…`, `92c997bc…`). The system's `gzip -t` and `zcat | md5sum` agree.

- Launcher: `sh /tmp/if2453/c27/bin/c27_fetch.sh launch B00 chain`.
- Each batch is its own unit, `c27-fetch-<kk>`: `systemd-run --user`, `nice -n 19`, `MemoryMax=8G`,
  `CPUQuota=200%`.
- At the end of a DONE batch the unit launches the next batch. Before every launch it checks STOP, that no other
  `c27-fetch-*` unit is active, and that `df /tmp` minus the batch's expected size (0.30 bytes per bp) stays at or
  above 250 GB.
- A batch that is refused, hits STOP, or has more than 5% of its genomes not served ends FAILED and stops the
  chain (B00 is exempt).

### Per batch (*measured*, `genomes/<B>/DONE`, `fetch_result.json`)

| Batch | Asked | Verified | Not served | Downloaded bytes | Kept `.fna.gz` bytes | FASTA bytes | Requests | Wall s |
|---|---|---|---|---|---|---|---|---|
| B00 | 7 | 0 | 7 | 24,171 | 0 | 0 | 8 | 4.1 |
| B01 | 2,406 | 2,396 | 10 | 2,887,732,717 | 2,885,360,389 | 9,971,590,448 | 131 | 597.5 |
| B02 | 2,639 | 2,639 | 0 | 2,897,518,626 | 2,894,797,085 | 10,029,249,403 | 132 | 608.1 |
| B03 | 2,381 | 2,381 | 0 | 2,885,443,272 | 2,883,179,085 | 10,024,685,232 | 120 | 567.4 |
| B04 | 2,357 | 2,356 | 1 | 2,896,262,694 | 2,893,986,120 | 10,018,070,116 | 119 | 583.2 |
| B05 | 2,376 | 2,376 | 0 | 2,907,795,831 | 2,905,558,586 | 10,021,035,912 | 120 | 647.0 |
| B06 | 2,240 | 2,240 | 0 | 2,912,684,634 | 2,910,576,756 | 10,025,030,041 | 112 | 591.7 |
| B07 | 2,365 | 2,365 | 0 | 2,905,973,477 | 2,903,772,408 | 10,025,219,391 | 119 | 568.3 |
| B08 | 2,286 | 2,286 | 0 | 2,910,269,964 | 2,908,202,564 | 10,024,415,755 | 115 | 539.9 |
| B09 | 2,341 | 2,339 | 2 | 2,907,196,759 | 2,905,037,609 | 10,019,228,143 | 120 | 567.6 |
| B10 | 2,290 | 2,290 | 0 | 2,910,224,866 | 2,908,110,355 | 10,028,033,691 | 115 | 545.7 |
| B11 | 2,278 | 2,274 | 4 | 2,895,558,113 | 2,893,555,132 | 10,002,163,928 | 118 | 523.1 |
| B12 | 2,524 | 2,515 | 9 | 2,898,818,702 | 2,896,484,985 | 9,976,252,796 | 136 | 494.9 |
| B01 to B12 | 28,483 | 28,457 | 26 | 34,815,479,655 | 34,788,621,074 | 120,164,974,856 | 1,457 | 6,834.4 |

- B01 had 0 refusals, 0 retries and 0 md5 mismatches, and every verified genome's bp equals GTDB's
  `genome_size`. The rates (*derived*) were 240 genomes per min and 4.83 MB/s downloaded.
- Kept bytes are 0.2930 per bp (*derived*, 2,885,360,389 / 9,847,934,765), close to phase 1's 0.29026 download
  ratio.
- The whole fetch (*measured*, the table's last row; sums *derived*): 34.79 GB kept against the 34.8 GB *estimate*; 113.9
  minutes of fetching (6,834.4 s), 21:23 to 23:17 CST; 118,675,230,482 bp verified, every genome's bp equal to GTDB's
  `genome_size` (0 mismatches in every batch); 0 md5 refetches.

**Not served.**

- B00, retried once, failed again: GCA_900637225.1, GCA_902153285.2, GCA_902502775.2, GCA_902502895.2,
  GCA_963667335.1, GCA_963921365.2, GCF_040719245.2. Datasets returned no zip for the batch, or for
  GCF_040719245.2 alone. They are left out (DEC-31).
- B01: 10 GenBank accessions, among them GCA_006228445.1, GCA_006228465.1, GCA_006228485.1, GCA_006228505.1,
  GCA_006228525.1, GCA_006228545.1, GCA_021559635.1, GCA_021559655.1 and GCA_021562195.1. The full list is in
  `genomes/B01/fetch_result.json` (`failed`).
- All 26 of B01 to B12 (*measured*; each "not returned by Datasets v2 (alone or in its batch)", the FTP fallback
  failing too): B01: GCA_006228445.1, GCA_006228465.1, GCA_006228485.1, GCA_006228505.1, GCA_006228525.1,
  GCA_006228545.1, GCA_021559635.1, GCA_021559655.1, GCA_021562195.1, GCA_040822005.1. B04: GCF_006228425.1. B09:
  GCF_028215935.1, GCF_028471605.1. B11: GCF_039544385.1, GCF_039545775.1, GCF_039602145.1, GCF_039614705.1. B12:
  GCF_041463605.1, GCF_041463615.1, GCF_041463625.1, GCF_041463635.1, GCF_041463645.1, GCF_041463755.1,
  GCF_041464225.1, GCF_041486685.1, GCF_041489455.1. Runs of consecutive accessions suggest whole submissions
  withdrawn. They are missing from the ANI candidates: 26 of 28,483 (0.09%).
- The likely cause, not checked, is suppressed assemblies.

**Disk** (*measured*, `logs/disk.log`):

- 560,492,540 KiB free on `/tmp` at the start.
- 556,177,660 KiB free at 21:34, with 3,205,898,656 bytes under `genomes/`.
- The peak is the kept total, about 35 GB (*estimate*), because nothing is unpacked. Other cards' builds share the
  disk, so the 250 GB floor is rechecked before every batch.

## State at handback (2026-11-12 21:34 CST) and how to resume

- Units: the chain runs one `c27-fetch-<kk>` at a time (02 was active; 03 to 12 follow by themselves).
- How to poll, with short ssh calls:
  - `tail -n 5 /tmp/if2453/c27/logs/launch.log`: `<B> exit <rc>` and `launch` lines; a `HALT` line means it
    stopped, with the reason, which is also in `genomes/HALT`.
  - `ls /tmp/if2453/c27/genomes/*/DONE /tmp/if2453/c27/genomes/*/FAILED`.
  - `systemctl --user list-units 'c27-fetch-*'`.
- Done: `genomes/B01` to `genomes/B12` each hold `DONE`, and no unit is active.
- To resume after a FAILED or HALT batch: `sh /tmp/if2453/c27/bin/c27_fetch.sh launch <that batch> chain`. A batch
  is resumable: genomes in its `lock.json` are not asked again, and DONE batches are skipped.
- To stop: `systemctl --user stop c27-fetch-<kk>`, or create the data root's STOP file.
- After it ended (done 2026-11-12 by the orchestrating session: the table, the `failed` list, `sources.yaml`); still to
  do, in step (4):
  - build `paths.tsv` with `c27_skani.py paths` from the phase-1 locks and B01 to B12's locks.

## Pod use

- Units: `c27-select` (1 CPU, 31 s) and the `c27-fetch-*` chain (`CPUQuota=200%`; about 0.6 CPU in use).
- Test runs: `pod_pytest.sh` under `/tmp/if2453/ci/C-27`.
- A few short scripts: split, selftest, dry run.
- Written only under `/tmp/if2453/c27` and `/tmp/if2453/ci/C-27`. The phase-1 locks under `/tmp/if2453/r50` and
  `/tmp/if2453/ladder` were only read.
- Network requests from the pod: only `c27-fetch-*`, to the listed Datasets URL. Load was about 48 (other cards).
