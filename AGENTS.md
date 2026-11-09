# Rules for agents working on bio-box

These rules bind every LLM agent that works in this repository. Tasks are in
[`planning/tasks/TASKS.md`](planning/tasks/TASKS.md). The evidence behind the rules is in
[`planning/measurements/2026-11-09_design_stress_test.md`](planning/measurements/2026-11-09_design_stress_test.md).

If a rule and a task card disagree, the rule wins and you stop and report the conflict.

## 1. How to work

1. Do exactly one task card per session. Read the whole card, then every file it lists, before editing.
2. Authority, highest first: `planning/DECISIONS.md`, the invariants in `00-spec.md` §3, this file, the task card, other plan documents. Plan documents contain superseded text; never act on a plan document against a decision row.
3. Work on a branch named `task/<ID>` (for example `task/B-03`). Never commit to `main`. Never push to the `pod` remote. Never force-push.
4. Edit only the files the card lists under **Files**. If the task cannot be finished without touching another file, stop and report which file and why.
5. Cards are tagged by tier. **L**: do it alone. **M**: do it, then request review before merge. **H**: not for you; a design or PI decision. If your card depends on an unfinished **H** card, stop.
6. Do not start a card whose **Depends** list has an unmerged card.
7. Do not redesign. If you think the card is wrong, finish nothing, and report the reason with evidence (a failing command and its output).
8. Shell: one command per call; no `&&`, `;`, `||` or pipes; quote paths that contain spaces; use `git -C <path>` outside the working directory.

## 2. Files you must not change unless the card names them

- `planning/DECISIONS.md`: append-only, and only an **H** card or the PI adds rows.
- `tests/registry/schema_lock.json`, `tests/registry/json_schema_lock.json`: never by hand; only by running `scripts/update_schema_lock.py` in a card that says so.
- `golden.yaml`, `golden.lock.yaml`, `refdb.yaml` checksums, `uv.lock`, `container/environment.yml`.
- `00-spec.md`, `PLAN.md`, `STAGES.md`, `DEPLOYMENT-AND-MALLEABILITY.md`.
- Any existing test, except to migrate it as the card describes. Never delete, skip, xfail or loosen a test to make it pass. Never add a retry.
- No new dependency unless the card names the package; then add it with `uv add` and nothing else.

## 3. Data model (hard rules)

1. Coordinates are 0-based half-open. Strand is `+1` or `-1`. On a circular contig a feature that crosses the origin has `end > contig_len`.
2. Every identifier comes from `src/biobox/ids.py`. UIDs are 16-byte BLAKE3 digests stored as `FIXED_LEN_BYTE_ARRAY(16)`; exact stored bytes are identified by sha256.
3. **An ordinal is never a key.** A position in a tape (`ordinal`, list index) is valid only inside one `tape_version`. It may not appear in a primary key, a foreign key, an exported file, or a search result without the `feature_uid` beside it. Features are keyed by `feature_uid = blake3(contig_uid, start, end, strand, kind)`. Evidence: with ordinal keys, adding one ncRNA to *E. coli* K-12 rewrote 1,214 of 4,318 slot rows and 24,342 of 86,360 pair rows; with coordinate keys it changed 0 and 110 (measured 2026-11-09).
4. Token ids are release-local. Never store a token id outside a `release=` partition. Never accept or return a bare token id at an API boundary; use the stable key (`model_acc`, `cluster_uid`) and resolve it through `token_dict`.
5. **Do not materialize pairs.** Neighbor questions are answered from `token_postings` (one row per slot and token) and the tapes. `pairs_of_contig` survives only as a test oracle. Evidence: `slot_pairs` measured 54 rows and 422 bytes per gene at three token levels, more than the sequence shards (334), and it grows faster than its budget (§3 of the measurements file).
6. Contig and genome attributes that can change without the sequence changing (taxonomy, tier, assembly quality, species-rep flag) live only in the dimension tables `contig_dim` and `genome_dim`. Never copy them into a fact table.
7. A table exists only as a `TableSpec` in `src/biobox/registry/schemas.py`. Files live only at the spec's `location`, partition keys appear only in the path, and every file is sorted by the spec's `sort_key` with min/max statistics on the sort-key columns.
8. Until the first release (`rS0.1`) is signed, a card may change a schema and regenerate the lock. After that, a locked column is never removed, renamed, retyped or reordered; new columns are appended and nullable (LY-09).
9. The manifest is the record. Readers use explicit file lists from the manifest, never globs. A file is written under a temporary name and renamed into place.
10. A stage never rewrites an earlier stage's file. This includes index and metadata tables: they are partitioned by stage, not rewritten.

## 4. Builders (hard rules)

1. **Bounded memory.** A builder may not hold a whole table, a whole stage's sequences, or a whole index in Python objects. It processes one genome or one batch at a time and writes through `biobox.io.write_table` (card B-01). Evidence: the pair builder used 544 bytes of heap per row and `SeqStore` held 414 MB and took 2 s to open per million contigs; S6 has 111.5 million contigs (measured 2026-11-09).
2. Every builder card has a memory test: the build of the 100× synthetic store (`scripts/stress/synth_scale.py`) stays under the peak RSS the card states.
3. **Deterministic.** The same inputs give byte-identical files on any host, at any thread count, in any input order. No timestamps, hostnames, absolute paths, random seeds from the clock, or dict-order dependence in an output. Sort inputs explicitly.
4. A lookup by key must prune. Any reader or SQL macro that resolves a key must read O(1) row groups, including when the keys come from a joined table. Evidence: the current sequence macros read 1 of 9 row groups for a literal key and 9 of 9 for a joined key (measured 2026-11-09). Test with `catalog.profile_query`.
5. Each step writes ledger rows with the versions of the tools that shaped its output.
6. Row groups: 122,880 rows for row tables; for list tables (tapes, tape tokens) a row group holds at most 131,072 slots in total, not a fixed number of contigs.

## 5. Tests (hard rules)

1. Write the test first and run it to see it fail. Put the failing output in your report. A test that never failed proves nothing.
2. An oracle is an independent implementation: an external tool (`samtools`, the Prodigal binary, `hmmsearch`, MMseqs2), a brute-force loop written separately, or a published table. It is never the function under test, and never a snapshot of that function's earlier output.
3. Fixed seeds. No sleeps to wait for state: poll a condition with a timeout.
4. Mark each test with `@pytest.mark.covers(...)` naming the plan test IDs it covers, and the ring marker.
5. Before you finish, run `uv run pytest -q` and put the summary line in your report. Tests that need MMseqs2, samtools or `hmmsearch` skip on a laptop; they run in the image in CI, and the card is not done until CI is green on your branch.
6. Statistics are tested only against constructed truth (planted signal recovered, null worlds calibrated). A golden snapshot is not evidence that a statistic is right.

## 6. Numbers

Every number you write in a document carries one of three labels: *measured* (host, date, command), *derived* (formula and its inputs) or *estimate*. Never copy a number from a plan document into code or into another document without the label and source. If a number you measure differs by more than 20% from the estimate a card or plan gives, stop and report it.

## 7. The pod (wilkcomp02)

It is a shared machine and we are guests. An agent runs something on the pod only when its card says **Pod: yes**.

1. Run only through `bio-box run-step` (card C-02), inside the `.sif` pulled by digest. Tools installed on the pod are never used to produce outputs; they may be used as oracles.
2. At most 40 threads by day; 56 at night and on weekends, and only when `ps` or `top` shows no other compute (the profile's 08:00 to 20:00 weekday window is an assumption until the PI confirms it with the pod owner). Every process runs inside a `systemd-run --user --scope -p MemoryMax=200G` scope, at `nice 19`. A single process stays under 300 CPU-minutes; larger work is split into shard processes.
3. Work in `/tmp/if2453/bio-box`. Results go to the data root with one `rsync` and an atomic rename. Never write files under 64 MB onto `/stor` except the registry's `small_ok` tables.
4. Before each step, check for a `STOP` file in the data root and exit cleanly if it exists.
5. Never delete or overwrite anything under `/stor`. Never use Docker on the pod. Never start a network service.
6. Code reaches the pod by `git push pod` from the PI's laptop. Agents do not edit files in `~/bio-box` on the pod.

## 8. Outside services

Fetch only URLs listed in `golden.yaml`, `refdb.yaml` or `sources.yaml`. Send no names, email addresses or project details in any request. Do not download or run a binary that is not pinned in `container/environment.yml` or `uv.lock`.

## 9. Stop and report instead of guessing

Stop when: the card is ambiguous; a listed file does not match what the card describes; a test outside your card fails; you would need to change a file from §2; a measurement misses an estimate by more than 20%; an acceptance command cannot be made to pass after two real attempts; anything on the pod behaves unexpectedly.

## 10. Report format

End every task with this report, and nothing else:

```
Task: <ID> <title>
Branch: task/<ID>   Commit: <sha>
Files changed: <list>
Test that failed first: <test id> and its failing output (first 10 lines)
Acceptance: <each acceptance command and its last line of output>
pytest summary: <line>
Numbers recorded: <number, label, source> or none
Deviations from the card: <none, or what and why>
Open questions: <none, or list>
```

Commit messages: `<type>(<scope>): <summary>` (types `feat`, `fix`, `docs`, `refactor`, `chore`, `test`, `perf`; imperative summary of at most 72 characters, no trailing period). Put the task ID and plan test IDs in the body. No footers, no sign-offs. Commit only; never push.
