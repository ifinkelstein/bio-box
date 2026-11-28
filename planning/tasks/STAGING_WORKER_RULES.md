# Standing rules for a staging-card worker (r26.11.5), from the orchestrating session, 2026-11-26

Read AGENTS.md in full first; it binds you. Then your card in `planning/tasks/TASKS.md`, the paragraph of common
rules above card C-73 with its dated line on DEC-80, your stage's record and sections 2, 5, 6, 7 and 8 of
`planning/measurements/2026-11-26_r26_11_5_staging_plan.md`, rows DEC-79 to DEC-83 of `planning/DECISIONS.md`, and
section W10 at the top of `planning/tasks/STATUS.md`. Read how the cards already done did the same things, and follow
them rather than invent: `planning/measurements/2026-11-26_mv2_mveukv_inputs.md` (C-79), `2026-11-26_egv1_inputs.md`
(C-69), `2026-11-26_sag2_inputs.md` (C-77), `2026-11-26_grare_inputs.md` (C-75), `2026-11-26_r232_index.md` (C-74),
with their scripts under `scripts/stress/build_stage_pod/`.

## Decisions already taken (do not reopen)
- Stage names have no hyphen. Confs are written as `stages/templates/<STAGE>.conf` only, with real `SUM=` lines.
  Nothing goes under `/tmp/if2453/s2/inputs/`. You do not edit `tests/r1/test_stage_driver.py`.
- DEC-80: a missing or unconfirmed licence does not hold a staging. Record the licence as an entry or record states
  it, or say none is stated; `redistributable: true` only where a recorded licence permits it, else `false`.
- No genome in hand is left out because GTDB R232 or GlobDB r232 lists it, or because r232 dropped it. Record its
  r232 status (the index's `lookup` or `rekey`: in GTDB R232 or not, representative or not, in GlobDB r232 or not,
  or "not visible to the accession test") in a table beside the stage's metadata, with counts by status in the note.
  Where your card says to remove on R232 or r232 membership, record instead. Removals are only: exact repeats within
  the set (DEC-27), what the genetic-code screen or `s2_check_inputs.py` refuses, non-sequence files, and what the
  card's own catch names for another reason (for example isolate reference genomes that are not the set's subject).
- An ANI dedup the card asks for runs against genomes that are on the pod: the staged stages and the GlobDB r226
  archaeal references of card C-75 (`/stor/scratch/Wilke/if2453/bio-box-staging/GRARE/globdb_r226_archaea.fna.gz`,
  list `/tmp/if2453/c75/meta/archaea/archaea_list.tsv`). No r232 genome is on the pod. An ANI match to a held stage
  is recorded (`EXPECT`/overlap tables); it removes a genome only where the card says the set is deduplicated against
  that held stage (for example `MIR1`'s GOEV archive against `GV1`).
- The index of card C-74 is `/tmp/if2453/c74/index/`. Only AMXMAG, GEM, GOMC and SPIRE carry INSDC accessions in
  GlobDB r232, so the accession test sees GTDB R232 in full and GlobDB r232 in part: report seen and unseen counts.

## Bytes
- Sets on the NAS reach the pod only with `scripts/orchestrator/relay_to_pod.sh` (`--dry-run` first; the directory's
  `SHA256SUMS` taken to the laptop with `scp` as the list; files of 64 MiB and more may go under
  `/stor/scratch/Wilke/if2453/bio-box-staging/<STAGE>/dl`, smaller ones under `/tmp/if2453/<card>/dl`). The script
  refuses a file name with a space: if a file you need has one, stop that part and report it; do not rename anything
  on the NAS and do not copy by hand.
- Read every archive's member list before staging and say in the note exactly what it holds. If it is not what the
  card expects, stage what is genome sequence of the card's subject, say so plainly, and leave the rest (with its
  reason) for the card it belongs to (eukaryote genomes: card C-60).
- The mirror of the staged inputs to `tower:/mnt/user/research-bulk/bio-box/sources/<STAGE>/` with
  `scripts/orchestrator/relay_to_tower.sh` comes last. Write nothing else on the NAS; never delete or overwrite there.

## The pod (`if2453@wilkcomp02.ccbb.utexas.edu`)
- It is shared, and a wave of 18 stage builds of ours runs on it (40 threads by day, up to 80 at night for some
  steps, under a guard that stops the build if a process passes 280 CPU-minutes). Units named `r26-11-4-wave-*`,
  `<stage>-<n>`, `<stage>-rh-*` and `<stage>-queue` are the build.
- Before each pod step look: `ssh if2453@wilkcomp02.ccbb.utexas.edu "systemctl --user list-units --type=service
  --state=active --no-legend --plain"`. If two staging units (names `c7*-*`, `c8*-*`) are active, wait a few minutes
  and look again; never start a third. Start no compute while a build unit of stage `MV1` or `V3` is calling genes.
- At most 8 threads in all for an ANI step and 4 for anything else, `nice -n 19`, in a detached unit
  (`systemd-run --user --unit=<card>-<what> --collect -p MemoryMax=32G -p MemorySwapMax=0 ...`); no process over
  300 CPU-minutes. Tools come from the image, as the earlier cards ran them.
- Read `/tmp/if2453/r10/root` only. New files of 64 MB and more go under
  `/stor/scratch/Wilke/if2453/bio-box-staging/<STAGE>/`; never overwrite or delete anything under `/stor`. Smaller
  files under `/tmp/if2453/<card>/`.
- Wait for a detached unit with one `systemctl --user is-active` call every few minutes. Do not end your work with a
  unit of yours still running.
- If something is missing, or a count differs by more than 20% from the card's, stop and report (AGENTS section 9).

## Tests
- Tests first (AGENTS section 5.1), on fixtures, seen failing, with the failing output in your report.
- Never run pytest or computation on the laptop. The only test command is
  `env BIOBOX_TREE=<worktree> sh <worktree>/scripts/pod_pytest.sh <card> <your test file> tests/r1/test_r232_index.py tests/r0`
  (one process; never the whole suite, never `--shards`).
- Tests must pass on GitHub Actions, which has no pod tools: a case that needs one skips cleanly. Do not write a test
  that reads `sources.yaml` directly (`tests/r0/test_image_files.py` rejects it). For the tools-hidden `PATH` do as
  `tests/r1/test_r232_index.py` does.

## Shell, files, commits
- One command per Bash call; no `&&`, `;`, `||` or pipes, in local commands and inside the quoted command of an ssh
  call alike (send a script file for a remote sequence). Quote paths with spaces. `git -C <path>` outside the
  working directory. If the session refuses a command, do not find another way to do the same thing: report it.
- Edit only the card's Files; in `sources.yaml` only the entries the card names. Every number carries a label
  (measured with host, date and command; derived with its formula; estimate).
- Do not touch plan/task-cards, `TASKS.md`, `STATUS.md` or `DECISIONS.md`. Commit on your branch only: message
  `<type>(<scope>): <summary>` (at most 72 characters), the card ID in the body, no footers, no sign-offs, no
  Co-Authored-By line. Do not push.

## Report
End with the report format of AGENTS.md section 10. Before it, at most six lines for the orchestrating session: per
stage, genomes, contigs, bp, bp new against held stages, counts by r232 status, and what was left out and why; where
the files are on the pod and on the NAS; pod compute used; anything that differed from the card.

## Ledger synchronization at fetch completion and merge (DEC-59, DEC-88, DEC-89)
- At every `tower` fetch completion and every staging-card merge, the orchestrating session (or the ledger card it
  assigns) synchronizes `sources.yaml` (`fetched:`, `stored:`, `fetched_pod:`, and `staged:`) and
  `planning/CATALOG_GOAL.md` (§2 order, §3 stage row matching every `.conf` in `stages/` and `stages/templates/`, §4.1
  `left_out.tsv` / deferred remainder, and §4.3 source status), and checks both files with
  `python scripts/orchestrator/r26_sync_sources_ledger.py verify` before closing the merge.
