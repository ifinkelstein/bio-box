# The complete genomes that S2C phase 2's selection left out, by accession (card C-51, DEC-59)

Date: 2026-11-21, host wilkcomp02 (the pod). Every number carries a label: *measured* (host, date, command),
*derived* (formula and inputs) or *estimate*. Card C-27 capped each GTDB R226 species at 20 complete genomes (DEC-34):
phase 1's genome and 19 picks, which stages `C1` to `C6` hold. Its candidate and remainder lists under
`/tmp/if2453/c27/` are gone (checked 2026-11-21), so this card makes the list again from the metadata the selection
read, and keeps it.

## Result

**23,866 complete genomes of the 243 species over the cap are left out** (*measured*, `c51_deferred_list.py`, 2026-11-21,
below). One row per assembly accession, sorted, with species and assembly level; 26 of them carry `not_served = y`.
It is [`2026-11-21_c51_deferred_list/left_out.tsv`](2026-11-21_c51_deferred_list/left_out.tsv) (1,374,419 bytes,
under the card's 5 MB, so committed whole):

- sha256 `d0f10fd5e7a0fdc6aa2b3513aa4caea8230ee459dc53d19f152ef4f2de4e55a1`, 23,866 rows and a header; columns
  `accession`, `species` (GTDB `s__` name), `assembly_level` (always `Complete Genome`: the rule keeps no other),
  `not_served`.
- Beside it: `species_check.tsv` (per species: counts, the table's figure, the stage, `differs`, `why`; sha256
  `d76227782aa1a5f8a83764e4aa86e2ca63fa7561385b0c4f025fe8af4d4e370b`), `summary.json` (counts, input digests) and
  `not_served_candidates.txt` (the 26, copied from `2026-11-13_s2c_skani_fetch.md`, "Not served").

**The counts agree with `ani_species.tsv` in every one of the 243 species (0 differ), and `C5`'s species leave 10,318
genomes, the figure card C-47 derived, exactly** (*measured*, `summary.json` `by_stage.C5` = 10,318; card C-47's
*derived* 11,358 − 988 − 52 = 10,318). No species is named because none differs: for each, the metadata's complete
genomes equal the table's `complete_genomes`, phase 1's genome is among them, and exactly 19 picks of `C1` to `C6` are.

## 1. The inputs

| Input | Where | sha256 |
|---|---|---|
| GTDB R226 `bac120_metadata_r226.tsv.gz` (235,663,508 bytes) | `/stor/scratch/Wilke/if2453/all-the-data/provisional/s1_inputs/gtdb_r226_metadata.tar`, read only, extracted to `/tmp/if2453/c51/gtdb/` | `45f1614aff80c19752860880b47b3cdeb66def8ff66dc106909b28bb2bd9bcbf` |
| GTDB R226 `ar53_metadata_r226.tsv.gz` (5,673,200 bytes) | the same tar | `f12310f223e6dd27dc942379de1cfbeec5ebd2b917786cf2ce28695555eab6e5` |
| the tar (241,346,560 bytes) | the same | `d686a3bf3541eeb4610798da39a0174c0661d15a3485ef19d850e10aab77a817` |
| phase 1's `all.tsv` (11,973 genomes) | `2026-11-11_r50_inputs/selection/all.tsv` | `e68d2576067d4d83f77a50384865de5470e6dff1a42716a06148cdf26c0efeaf` |
| `stages.tsv` (13,477 picks of `C1` to `C6`) | `2026-11-15_s2c_phase2_inputs/stages/stages.tsv` | `d15649b2b5157a00c68204f971e3119af1b2621b00a0d1b9ac11cb02025cb49d` |

*Measured* (`sha256sum` on the pod, 2026-11-21): the two metadata files and the tar equal the values in `sources.yaml`
(`gtdb_metadata`), and the two files equal the `gtdb_files_sha256` that card C-27's `selection.json` recorded when it
read them, so this is the table the selection read. It was **not fetched again**: `sources.yaml` gives only GTDB's
release directory, not a file URL, and says no fetch is approved until an exact URL is added; the copy on `/stor` made
the fetch unnecessary. No network request was made. `/stor` was read, never written.

## 2. The rule and the script

`scripts/stress/c51_deferred_list.py` (stdlib, `python -I`, in the image `bio-box-test_faec1f4.sif`, `nice 19`, one
process, a few seconds): the complete genomes (`ncbi_assembly_level = 'Complete Genome'`, as `c25_select.py`) of every
species with more than 20; remove phase 1's genomes and the picks of `C1` to `C6`; write the rest. It refuses a
not-served candidate that is also a pick, a species in the table that is not over the cap in the metadata, and a
species with picks in two stages. Output has no timestamp or path and does not depend on input order.

Test first (`tests/r1/test_c51_deferred_list.py`, a world of 6 species at cap 2 worked by hand: a species of 5 complete
genomes and a scaffold, one at exactly the cap, one with a not-served candidate, one whose phase-1 genome is not
complete in the metadata, one in another stage, one with a failed pick). With the script absent, all four tests
failed (`can't open file '.../scripts/stress/c51_deferred_list.py': [Errno 2] No such file or directory`); after it,
`4 passed` on the pod (*measured*, `pod_pytest.sh C-51`, 2026-11-21).

Run (*measured*, wilkcomp02, 2026-11-21 02:4x CST, one scope with `MemoryMax=8G`):
`python -I c51_deferred_list.py /tmp/if2453/c51/gtdb all.tsv stages.tsv /tmp/if2453/c51/out --cap 20 --species-table
ani_species.tsv --not-served-candidates not_served_candidates.txt --not-served-picks not_served.txt`.
The metadata has 49,324 complete genomes in 11,973 species; 243 species have more than 20, with 28,726 complete
genomes (card C-27's candidates: 28,726, equal).

## 3. The check against `ani_species.tsv`

- *Derived* total: 28,726 complete genomes − 243 × 20 kept = 23,866 left. *Measured*: 23,866 rows. Equal.
- Cross-check against card C-27's own figures: the fetch kept 28,457 of 28,483 candidates (23,840 not kept, removed
  from `/tmp` 2026-11-14) and 26 were not served: 23,840 + 26 = 23,866 (*derived*). Equal.
- Per species (`species_check.tsv`): `left` equals `expected_left` (the table's `complete_genomes` − 20) in all 243
  rows; `phase1_in` is 1 and `picks_in` is 19 in all 243. `species_differ` = 0, `picks_not_complete_in_metadata` = [].

Per stage, by the stage that holds the species' picks (`stages.tsv`; a species is whole in one stage, asserted):

| Stage | Species over the cap | Genomes left out | of which not served (candidates) | Picks not served in its species (the 8) |
|---|---|---|---|---|
| `C1` | 19 | 1,018 | 0 | 0 |
| `C2` | 67 | 5,979 | 2 | 1 |
| `C3` | 34 | 1,712 | 0 | 2 |
| `C4` | 33 | 2,116 | 0 | 0 |
| `C5` | 52 | 10,318 | 7 | 1 |
| `C6` | 38 | 2,723 | 17 | 1 |
| none | 0 | 0 | 0 | 3 |
| all | 243 | 23,866 | 26 | 8 |

(*measured*, `summary.json`; the totals are sums, *derived*.) For `C5`: 52 species, 11,358 complete genomes in them
(*derived*: 10,318 + 988 + 52, where 988 are `C5`'s picks and 52 are phase 1's), which is card C-47's figure to the
genome, so the assumption that card made ("every one of the 52 species' phase-1 genomes is among its complete genomes")
holds, checked genome by genome here.

## 4. The genomes NCBI did not serve (DEC-34 (3): 34 accessions)

Two different things, handled differently:

- **The 26 candidates** (card C-27, `2026-11-13_s2c_skani_fetch.md`): they were left out of the ANI pick for want of a
  sequence, so they are among the 23,866, **kept in the list and marked `not_served = y`**. By the stage of their
  species: `C2` 2 (*Streptococcus dysgalactiae*), `C5` 7 (*E. coli* 4, *Klebsiella variicola*, *Vibrio cholerae*, *V.
  harveyi*: one each), `C6` 17 (*Xanthomonas euvesicatoria* 8, *X. cissicola* 7, *X. arboricola* 1, *Acinetobacter
  nosocomialis* 1); the species are card C-27's nine, and which stage holds each is *derived* here from `stages.tsv`.
  Fetch them again when NCBI serves them.
- **The 8 picks that failed** (`2026-11-13_s2c_phase2_fetch/not_served.txt`): they were picks, in species at or under
  the cap, never fetched and in no stage. They are **not in the list** (a species at or under the cap has no remainder
  here). By the stage of their species: `C2` 1 (*Clostridium butyricum*), `C3` 2 (*Lactobacillus* sp945980025,
  *Neorhizobium banfieldiae*), `C5` 1 (*Klebsiella planticola*), `C6` 1 (*Xanthomonas hydrangeae*), and 3 in species
  that have no stage because the failed genome was their only remainder (*Xanthomonas codiaei*, *Neorhizobium
  flavum*, *Thiodubiliella endoseptemdiera*). Their register row stays as it is, with these counts added.
- Not covered here: the 7 phase-1 genomes NCBI did not serve (`2026-11-11_r50_inputs.md`) are phase 1's own, and phase
  1 is removed from the list, so they are in neither.

## 5. The copy on `tower`

The list is on the pod at `/tmp/if2453/c51/relay/` for the relay: `gtdb_r226_complete_left_out.tsv` (1,374,419 bytes,
the list above) and `gtdb_r226_complete_left_out_species_check.tsv` (36,394 bytes). The list file for
`scripts/orchestrator/relay_to_tower.sh` (`sha256  name`):

```
d0f10fd5e7a0fdc6aa2b3513aa4caea8230ee459dc53d19f152ef4f2de4e55a1  gtdb_r226_complete_left_out.tsv
d76227782aa1a5f8a83764e4aa86e2ca63fa7561385b0c4f025fe8af4d4e370b  gtdb_r226_complete_left_out_species_check.tsv
```

`relay_to_tower.sh --dry-run /tmp/if2453/c51/relay /mnt/user/research-bulk/bio-box/deferred_lists/C <that file>` (laptop,
`BIOBOX_POD_JUMP` set, 2026-11-21) printed `copy` for both files and `dry run: 2 files to copy, 1410813 bytes, about 0.0
hours ... 0 in place; 0 missing on the pod`. **The transfer was not run**: the orchestrator runs it after the audit.
The commit holds the same list, so the list is also in the repository; the copy on `tower` is the second one.

## 6. Open questions

- The file URL of GTDB's metadata is not in `sources.yaml`. If the list must be made again from GTDB itself, the PI
  or an `H` card adds it; until then the `/stor` tar (sha256 above) is the only copy of the table.
- The list names the genomes, not their sequence: a card that fetches them goes by accession from NCBI datasets as
  card C-27 did, and stages them after `C6` with `contig_uid` dedup (complete genomes carry plasmids that `PL1`
  stores).
