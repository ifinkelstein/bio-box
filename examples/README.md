# Worked examples

A worked example is a real question asked of a release with the commands of `GUIDE.md`, kept so that anyone can run
it again: one directory per example, one script that runs it on the pod, and the facts its answer must show. Before
a milestone release is announced, every example is run on it (AGENTS §5.7, card E-11; `RUNBOOK.md` section 6). A
failed invariant stops the announcement.

| Example | The question | Anchors | Cost on `rN.3` |
|---|---|---|---|
| [`rnap/`](rnap/README.md) | Where are the T7-like RNA polymerases, in phages and prokaryotes, and what sits next to them? (`GUIDE.md` sections 3 to 7, as a check) | Pfam `RNA_pol` (PF00940) at 30 bits; control T7 gene 1 | 18.2 and 18.4 CPU-minutes at 2 threads in two runs (*measured*, wilkcomp02, 2026-11-13, its README section 4 and `planning/measurements/2026-11-13_release_examples_rN3.md`) |
| [`crispr_cas_neighbourhoods/`](crispr_cas_neighbourhoods/README.md) | Minimal type I-F CRISPR-Cas systems carried by Tn7-like transposons (Peters et al. 2017) | Pfam `TniQ` and `Cas_Csy3` at their gathering cutoffs | 22.8 and 22.85 CPU-minutes at 2 threads in two runs (*measured*, wilkcomp02, 2026-11-13, its README section 6 and `planning/measurements/2026-11-13_release_examples_rN3.md`) |

## What an example holds

| File | What it is |
|---|---|
| `README.md` | the question, the plan written before the run, how to run it, what it found on a release, what each failed invariant would mean, its cost |
| `run.sh` | the whole example as one command on the pod host, in the environment below |
| `expected.tsv` | the checks: `check_id`, `asserts` (in words), `kind`, `expectation` |
| `recorded/<release>.tsv` | the recorded values of one release (`check_id`, `value`), copied from a run's `recorded_new.tsv`, never typed |
| `sql/` | the SQL the script runs (`bio-box sql`), and any SQL a table of the README came from |
| `results/` | the recorded run's small files: `check.txt`, `observed.tsv` (every value the check read) and the like |

Shared by every example, in `examples/_shared/`:

| File | What it is |
|---|---|
| `launch.sh` | sourced by each `run.sh` and by `scripts/release_examples.sh`: `room` waits until the account's CPU (the sum of `ps -u $USER -o pcpu=`) plus the launch's threads fits under the hour's cap (40 by day, 56 at night and weekends; AGENTS §7.2) and stops at a `STOP` file in the root; `bb NAME ARGS...` runs one `bio-box` command in its own capped unit through `scripts/pod_run.sh` and waits for it; `tool LOG ARGS...` runs one command of the image (HMMER, Python) in a capped scope |
| `check.py` | compares `observed.tsv` with `expected.tsv` and the release's `recorded/` file (standard library only); exits 1 when an invariant fails or is missing |
| `oracle.py` | the `hmmsearch` oracle: a run's hits against HMMER run by hand on chosen proteins, at a cutoff in bits or the model's GA |

They are tested by `tests/r1/test_examples_check.py` and `tests/r1/test_release_examples.py`; every `bio-box` command
of each `run.sh` and of each README's command blocks is parsed by `tests/r0/test_runbook.py`, so a renamed flag fails
a test before it fails a release check.

## The two kinds of check

- `invariant`: a fact that must hold on every release that holds the example's control genomes, whatever else the
  release holds: the control is found at or above the cutoff with its architecture, its known neighbours sit at their
  offsets, the search equals `hmmsearch` run by hand, the negative control is empty. Written as `== V` (equal as text,
  or as numbers when both parse), `>= N`, `<= N` or `in A,B,...`. A failed or missing invariant fails the check.
- `recorded`: a value of one release (a count, a name, a checksum). Its expectation is the word `recorded`; the value
  is compared with `recorded/<release>.tsv`. A difference is reported (`differs`), not failed: a new release is
  expected to differ. A rerun of the same release that `differs` means the release or the code changed: find out which
  before recording anything.

Among the invariants, each example has at least one **guard**: a check that the window sees what it should (genes known
to be there are there), so that an invariant of the form "none of X near Y" cannot pass because the window is empty.
And at least one check by SQL that does not go through `nbr` (coordinates and the release's Pfam scan alone), so that
a failure says whether the data or the tool changed.

`check.py` prints one line per check, `ok`, `FAIL` (an invariant), `differs` (a recorded value), `new` (no recorded
value for this release yet), `missing` or `unexpected`, then a summary line.

## Running one example

On the pod host, in `bash`, as `GUIDE.md` section 2 sets it up. Every example takes the same environment: `TREE` (a
checkout of current code, which the image runs as its `PYTHONPATH`), `PROFILE` (the host profile), `IMAGE` (the
`.sif`; it must hold matplotlib, which the runtime image built from `main` does), `DIGEST` (its sha256), `ROOT` (the
data root, bound read-only), `RELEASE`, and `OUT` (a new or empty scratch directory outside the root). Optional:
`THREADS` (2), `CAP`, `BIND`, `UNIT`, `MEM`, `PFAM` (each `run.sh` says what they do).

```bash
export TREE=/tmp/if2453/e11/tree
export PROFILE="$TREE/profiles/wilkcomp02.yaml"
export IMAGE=/tmp/if2453/images/bio-box_64f1bf0.sif
export DIGEST=sha256:cb7353f85316981be5268c2212100322b5fbf1086c32b07e2f9b47802b6afba3
export ROOT=/tmp/if2453/r10/root
export RELEASE=rN.3
export OUT=/tmp/if2453/<you>/rnap-rN.3
bash "$TREE/examples/rnap/run.sh"
tail -n 1 "$OUT/check.txt"
```

`run.sh` exits with `check.py`'s code: 0 when every invariant holds, 1 when one fails or is missing, 2 on a malformed
input or a step that failed (its log is named).

## Running all of them on a release

`scripts/release_examples.sh RELEASE` runs every example in turn, read-only on the root, each launch inside the
thread cap, into `$OUT/<example>/`, then checks each one and writes one summary, `$OUT/summary.tsv` (example, check,
kind, expected, observed, verdict), and a measurement-note skeleton, `$OUT/note.md`, to finish by hand into
`planning/measurements/`. It exits 1 when any invariant of any example fails or is missing. `--dry-run` prints what it
would run; `--check-only` runs only the checking step on an `$OUT` that already holds each example's `observed.tsv`.
So that it outlives the shell, start it as a small user unit of its own (each step is its own capped unit):

```bash
export OUT=/tmp/if2453/<you>/examples-rN.3
bash "$TREE/scripts/release_examples.sh" --dry-run rN.3
systemd-run --user --unit ex-driver --collect --quiet -p MemoryMax=2G -p MemorySwapMax=0 -p StandardOutput=append:/tmp/if2453/<you>/examples-driver.log -p StandardError=append:/tmp/if2453/<you>/examples-driver.log --setenv=TREE="$TREE" --setenv=PROFILE="$PROFILE" --setenv=IMAGE="$IMAGE" --setenv=DIGEST="$DIGEST" --setenv=ROOT="$ROOT" --setenv=OUT="$OUT" --setenv=USER="$USER" /usr/bin/nice -n 19 bash "$TREE/scripts/release_examples.sh" rN.3
tail -n 1 /tmp/if2453/<you>/examples-driver.log
```

Then: an invariant that fails is a finding to report, with the example's table of what a failure would mean (each
README has one); a `differs` on a new release is expected, and the release's own file is added by copying each
`$OUT/<example>/recorded_new.tsv` to `examples/<example>/recorded/<release>.tsv` in the commit that records the run,
with the finished note. Never edit another release's file.

## Adding an example

1. Choose a question whose answer is known independently (a paper, a genome map) and that the release can hold. Write
   the plan in the README before running anything: the anchors and cutoffs, the controls (a gene by
   `ACCESSION:START-END`, its `feature_uid` and `protein_uid` from `bio-box find`), what would count as agreement,
   and what bio-box cannot do for the question today.
2. Write `run.sh` in the form of the two here: the environment above, `. "$TREE/examples/_shared/launch.sh"`, each
   `bio-box` command as `bb NAME ...`, each other tool as `tool LOG ...`, nothing written under `ROOT`, then
   `observed.tsv` and the shared `check.py` with `--write-recorded "$OUT/recorded_new.tsv"`.
3. Write `expected.tsv`: invariants that hold on any release with the controls (with a guard, a check without `nbr`,
   and the `hmmsearch` oracle through `_shared/oracle.py`), and recorded values. The check ids must be unique; each
   observed value comes from SQL, the oracle or a script that reads the run's files, never typed.
4. Run it on the current release; copy its `recorded_new.tsv` to `recorded/<release>.tsv`, its `check.txt` and
   `observed.tsv` to `results/`, and write what it found in the README, every number labelled (AGENTS §6).
5. Add its invariant names to `tests/r1/test_examples_check.py` in the way the two examples' are; `tests/r0` parses
   its `run.sh` and README by itself. `scripts/release_examples.sh` finds it by its `run.sh` and `expected.tsv`.

## The rule (AGENTS §5.7)

"**Worked examples are release checks** (PI, 2026-11-13). Before a release is called a milestone (announced in
`STATUS.md`, `README.md` or `GUIDE.md` as the one to use, or copied to durable storage as such), run every worked
example under `examples/` on it, as a user would, read-only [...]. Each example's `invariant` checks (its positive
controls and structural facts) must hold; its counts are recorded per release in the example's `recorded/` directory
and in a measurement note, never edited for an earlier release. A failed invariant stops the announcement: report it
with the command and its output. An example is a check of the whole path from search to neighbourhoods on real data,
not a substitute for the tests of this section." (`AGENTS.md` §5.7.)

In practice: `scripts/release_examples.sh RELEASE` is the one command; `$OUT/summary.tsv` and its exit code are the
result; `$OUT/note.md`, finished by hand, is the measurement note; the `recorded_new.tsv` files become the release's
`recorded/` files.
