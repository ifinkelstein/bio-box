"""`bio-box build-stage`: a stage and its release as a fixed list of `run-step` steps (card C-04).

    bio-box build-stage --profile P --root R --stage S --release REL [--stages S0,S1]
        --genomes genomes.json --objects objects.json [--library NAME VERSION HMM CUTOFF ...]
        [--tmp D] [--threads N] [--scan-jobs J] [--scan-cpus C] [--cluster-parts N] [--gene-call-version V]
        [--ncrna-cm Rfam.cm.gz [--ncrna-set prok-1]]

No workflow engine: the list `steps` returns and the run records `run-step` writes. Each step is
`runstep.run_step` running `python -m biobox.cli stage ...` or `python -m biobox.cli release ...`
with `--out` left to `$BIOBOX_STAGING` (card C-08), in this order:

    inputs, ingest, genecall, features, igr, tapes, domains-<lib>-<ver> (one per --library),
    cluster, tokens, postings, dim, manifest, catalog

`features` is the slot for feature producers other than the gene caller. Without `--ncrna-cm` it has no
command and is passed over. With `--ncrna-cm F` (card E-07) it is the ncRNA producer
(`biobox.l1.ncrna`, model set `--ncrna-set`, default `prok-1`), a step named `ncrna-<producer>`, and the
IGR, tape, token and posting steps are given both producers (`--producers <gene caller>,<ncRNA
producer>`), which makes a tape version of its own by design; the gene-call files are the same with or
without it. The ncRNA scan runs with `--cpus` C (`--scan-cpus`, default `DEFAULT_SCAN_CPUS`, at most
`--threads`). A domain or ncRNA scan over the CPU budget runs as several steps (below). The stage
steps are recorded under `runs/stage=<S>/`, the release steps under `runs/stage=<REL>/`.

Contigs shared across stages (card C-20, DEC-26). Ingest stores a sequence once, in the first stage that sees
it, and that stage owns the contig's genes, IGRs and tape. When gene calling starts, `build-stage` adds
`--seq-stages` naming the stage and every other stage of the root (`with_seq_stages`), so a genome that holds
a contig stored by an earlier stage is still read and trained whole; gene calling skips that contig and
records the skip. Before the release's clustering (after its CPU estimate, below), `dim.owner_problem` refuses
(exit 2) a release that lists a stage without the stage that owns one of its contigs, naming that stage.

The step list and every argument are a pure function of the command line and the profile (the gene-call
step's `--seq-stages`, added when it starts, is the one exception). Paths
are made absolute. The thread count is `--threads`, else the smaller of the profile's day and
night caps, which `run-step` accepts at any hour; it is given to `run-step` (`BIOBOX_THREADS`), to
gene calling (`--workers`), to the scans (`--cpus`) and to MMseqs2 (`--threads`). Scratch is
`--tmp`, else `<profile tmp>/scratch`. Each step that takes a `--tmp` gets its own directory there,
`<scratch>/<stage or release>.<step>` (card C-13): before each attempt of the step (no run record, or
an intent file) build-stage removes what an earlier attempt left in it, and once the step exits 0 it
removes the directory. A killed builder never cleans its scratch (C-10: 1.2 GB after one stopped
clustering step, *measured*). The scan shards, their rescans and the merge share `<scratch>` itself:
the shards' runs there are the merge's input (below). The release steps take the release's stages from its
`release_def`: `--stages` is passed to the builder release steps only when it is given. On first
use it defines the release (the cluster step writes the definition), and later the steps verify it.

Before any step (exit 2, nothing run): a bad name; `--stage` not among the release's stages; no
`--stages` and no `release_def`; `--stages` that disagree with the `release_def`; a source object
whose own ingest estimate exceeds the profile's `cpu_min_per_process` (one process reads each object
whole, also in batches, card B-30); and `--scan-jobs` × `--scan-cpus`, or `--threads`, above the threads `run-step` allows at
the start of the call (`runstep.choose_threads`: the profile's cap for the hour and weekday, the day
cap while the host is busy). The ingest estimate is skipped when ingest's run record exists. Then
clustering's estimate, when every stage of the release has a gene-call run record and the cluster step
none (below). Before each step: a `STOP` file in the root exits 75, as `run-step` does.

Clustering's budget (cards B-36, B-42; DEC-29). The release's clustering in one process is estimated from the
rows of the release's `protein_seq` files (`cluster_cpu_min`, `release_proteins`): before any step when every
stage is gene-called, else (a stage of this call is not yet) just before the cluster step's command. Above
`SHARD_MARGIN` of the budget the release clusters in parts (below), so in one process it is within the budget;
should it ever be over `cpu_min_per_process` there, it is refused (exit 2), naming DEC-29 and the way to cluster in
parts, whatever `cluster_cpu_min_max` says. DEC-25's exemption of the one-process clustering step ended with card
B-42 (clustering in parts reproduced `rN.1b` byte for byte, planning/measurements/2026-11-12_split_clustering_rN1b.md);
the one exemption left is DEC-29's, for a call that no number of parts fits (below), under the profile's
`cluster_cpu_min_max`. Every other step keeps the 300 CPU-minute rule. A `STOP` that appears while a step
runs lets the step's command finish; `run-step` then publishes nothing of it and exits 75, so the next
call runs that step again from the start (the documented rule of card C-02; its cost in C-10: S1a's
gene calling, 57.7 s and 4.98 CPU-minutes, ran twice, *measured*).

The hour (card C-12). A build can run past the start of the day window, so the allowance is not
decided once: it is read again (`clock`, `load`, `runstep.choose_threads`) before each step that has
no run record and before each scan shard launch. Work is started only while it fits: a step only
while its `--threads` is at most the allowance, a shard only while the shards in flight and it, at
`--scan-cpus` each, are at most the allowance. Work already running finishes. When the allowance has
shrunk below `--scan-jobs` × `--scan-cpus`, fewer shards run at once; below one shard's
`--scan-cpus` (or a step's `--threads`), the launch is held: build-stage waits on the shards it
runs, or for `RECHECK_S` seconds when none runs, then checks again (and checks `STOP`). The default
`--threads` (the smaller cap) is never held. `run-step` itself chooses its threads by the hour when
it starts. The output bytes depend on none of this.

A scan in shards (cards B-24, C-11). Residues are known only after gene calling, so a scan step's
estimate is made when the build reaches it and its record does not exist. When it exceeds
`SHARD_MARGIN` of the budget (or a shard plan for it is already published), the scan runs as steps
of its own: `domains-<lib>-<ver>.plan` (`bio-box stage domains plan`, shards of at most
`shard_residues(profile)` residues: what fits `SHARD_MARGIN` of the budget at `SCAN_CPU_MIN_PER_M`),
then `domains-<lib>-<ver>.shard-<k>` for each shard (`... shard --shard K --cpus C`), then the merge
(`... merge`) under the scan step's own name, so its run record is the scan's. The plan depends on
the inputs, the budget and the margin only, never on J or C, and the merged files are those of the
scan in one process, byte for byte. The shards' scratch runs live under the scratch directory and
are removed after the merge.

Ingest in batches (card B-30). When ingest's estimate in one process exceeds `SHARD_MARGIN` of the
budget (or a pass-1 record of the stage is already published), ingest runs as steps of its own
(`ingest_steps`, `biobox.l0seq.store`): `ingest.pass1-<k>` for each pass-1 shard (the objects by name,
cut at `pass1_bytes`: what fits the margin at `INGEST_CPU_MIN`), `ingest.plan`, `ingest.batch-<j>` for
each batch of the plan (one shard each, DEC-21), then the merge under the step name `ingest`. Each is a
one-thread `run-step`; the pass-1 and batch steps run up to J × C at once within the allowance, as scan
shards do. The files are those of ingest in one process, byte for byte; the pass-1 and plan records
(`l0seq/ingest_plan/`) are a build's own, never in a manifest. Each of these steps has its own scratch
directory (`--tmp`, cleaned as above); the pass-1 runs and the decision table they hand on live in
`<scratch>` itself (`--runs-tmp`), are made again when lost, before the plan and before the batches
(`_rescan_lost_ingest`, as below), and are removed after the merge.

Clustering in parts (card B-41, `biobox.l2.cluster`). When the release's clustering in one process is
estimated over `SHARD_MARGIN` of the budget (`cluster_cpu_min`), or `--cluster-parts N` forces it, or a
plan of it is already published, the cluster step runs as steps of its own (`_cluster_in_parts`):
`cluster.plan` (the FASTA and `createdb`), then one step per module call of the fixed list
(`l2cluster.calls()`), a per-query call as `<call>.cut` and, when cut into more than one part,
`<call>.part-<p>` for each part (run up to J at once with C threads each within the allowance, as scan
shards) and the join under the call's own name; then `cluster` (`release cluster finish`: the tables)
under the cluster step's own name, so its run record is the clustering's. Parts are cut at
`SHARD_MARGIN` of the budget from the call's measured share (or into N), never from J, C or the hour;
the files are those of clustering in one process, byte for byte; the records (`l2/cluster_plan/`) are
a build's own, never in a manifest. The MMseqs2 databases live in `<scratch>/biobox-cluster/`; a join's
parts are removed once it is published, the work directory once the tables are. A call that finds the
plan published first makes again every published step whose outputs are not in scratch, in order
(`_redo_lost_cluster`, `bio-box release cluster redo`), each of which must reproduce its record.

One process for a call that cannot be split (card B-45, DEC-29). A prefilter whose parts each rebuild the
target's k-mer table (no shared index with `--diag-score 0`, card B-44) cannot be cut once that table alone
fills the part budget. With the profile's `cluster_cpu_min_max`, which the plan step is given as
`--one-process-cpu-min-max` and records in its cut, such a call's cut makes it one part, the workflow's own
call, and refuses it when its estimate (work plus table) is over that cap; without the key the cut refuses it.
Before the call's step is launched its estimate is checked against the profile again (`one_process_refusal`,
exit 2 above the cap or without the key), and the step runs with `budget_reason` `exempt: DEC-29`. The
release-level check above still comes first: in one process the release's estimate against the budget; in parts
(cards B-45, B-42), with or without the key, each call the published cuts run as one process, as before its launch
(`one_process_refusal`), and never the release's estimate in one process, so a release over the budget whose calls
all split under it is not refused. A call cut later is checked by its cut and again before its launch.

Lost scratch (card C-12). A reboot clears `/tmp`; a shard's run record then still exists while its
run in scratch is gone, so `run-step` would never scan it again and the merge would refuse. Before
the merge, every shard's run is checked against its published record (`domains.run_problem`: missing,
another sha256, or a record without one), and exactly those shards are scanned again into scratch
(`python -m biobox.l1.domains rescan`, each its own process under the cap and nice of a run-step,
run like the shards: up to J at once within the allowance, `STOP` before each launch). The records
stay as published: a rescan must reproduce its record's sha256, and one that does not stops the
build before the merge with an error naming the shard.

The ncRNA scan in shards (card E-07). Its estimate is the stage's genome bases (the `contig_len` of the
default variants in the stage's `contig_index`, listed in ingest's run record) at `NCRNA_CPU_H_PER_10GBP`;
over `SHARD_MARGIN` of the budget (or with its plan published) it runs as `ncrna-<p>.plan` (`bio-box stage
ncrna plan`, shards of at most `ncrna_shard_bases(profile)` genome bases), `ncrna-<p>.shard-<k>` and the merge
under the step's own name, exactly as a domain scan does, with the same pool, hours, STOP and lost-scratch
rescans (`python -m biobox.l1.ncrna rescan`). The pool has one more term, memory: each `cmsearch` holds up to
`ncrna_job_bytes(C)`, so at most `cap // ncrna_job_bytes(C)` shards run at once under the memory cap this
process runs in (the `pod_run.sh` unit's, `runstep.detect_cap`); a cap below one shard is refused (exit 2)
before the scan starts. The shards and the one-process step are given `--seq-stages` (every stage on disk) when
they start, as gene calling is, so a contig another stage owns is skipped and recorded (DEC-26).

Up to J shards run at once (`--scan-jobs`, default: the allowance at the start // C; fewer while
the allowance of the hour is smaller, above), each its own `bio-box
run-step` process (`python -m biobox.cli run-step --threads C -- ...`, `--scan-cpus`, default
`DEFAULT_SCAN_CPUS` or the allowance if smaller), started in this process's group, so a Ctrl-C or the
unit's stop reaches it and run-step forwards it to its step. Each launch is a thread that waits on
its process; build-stage waits on those threads and never polls. In shard order: a shard whose run
record exists (and no intent file) is passed over; before each launch a `STOP` file stops launches;
a shard that exits non-zero stops launches. Shards already running finish and publish either way,
then the call returns the failed shard's code (of the lowest shard when several fail), else 75 for
`STOP`. The next call runs only the shards without a run record, then the merge.

Concurrent run-steps into one root are safe because each has its own step name, so its own staging
directory, `.incoming/<stage>.<step>/`, intent file and run record; the shard records are distinct
files of a `small_ok` artifact (the size check does not group them) and the shards write no ledger.
The one directory they share, `<root>/.incoming`, `run-step` itself creates and removes safely beside
other steps (card C-12, `runstep._remove_leftovers`); build-stage keeps nothing there.

The first step that does not exit 0 ends the build with its exit code. `run-step` gives every step
kill-and-resume (cards C-02, C-03), so the next call skips each step whose run record and outputs
exist and resumes at the step that stopped. A second call after a complete build starts no
command. `build-stage` runs inside the memory-capped unit that `scripts/pod_run.sh` opens on the
host, and every step's `run-step` inherits that cap. It never calls `pod_run.sh` and never polls a
unit.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import json
import math
import os
import queue
import shutil
import subprocess
import sys
import threading
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq
import yaml

from biobox import catalog, dim, runstep
from biobox.l0seq import store
from biobox.l1 import domains, ncrna
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l2 import cluster as l2cluster
from biobox.profile import Profile

# Measured, wilkcomp02, 2026-11-10: one ingest process took 42.9 CPU-minutes for 9,596,180,328 bp of
# FASTA (stage P2 of the import ladder, planning/measurements/2026-11-10_import_ladder/README.md,
# "Step timings" and finding 4). The estimate applies the rate to bytes of decompressed FASTA, which
# are at least the bp they hold (headers and newlines add to them). The estimate therefore errs high.
INGEST_CPU_MIN = 42.9
INGEST_BP = 9_596_180_328
# Measured, wilkcomp02, 2026-11-10 (card C-10, planning/measurements/2026-11-10_build_stage_pod.md): the
# fixed start-up cost of a step, 0.07 CPU-minutes for stage S1a's `inputs` step, which does almost no
# work. Added to every ingest estimate (card C-13). Without it the estimate of S1a's ingest (63,357,113
# bytes) was 0.28 CPU-minutes against 0.36 measured; with it 0.35 (derived), and P1's (1,029,384,352
# bytes) 4.67 against 4.44 measured (derived). The rate above was measured with this cost included, so at
# P2's size the estimate is 0.07 high (0.2%, derived).
INGEST_FIXED_CPU_MIN = 0.07
# Measured, wilkcomp02, 2026-11-10: the Pfam scan of stage P1 in the first build with build-stage (card
# C-10, same file: the stage's 918,765 proteins, 294,516,869 residues, in 19 shard processes) took
# 2,343.74 CPU-minutes, 7.96 per million residues (derived: 2,343.74 / 294.52). The import ladder's scan
# of the 918,525 proteins new in P1 had measured 7.30 (2,149.9 CPU-minutes, planning/measurements/
# 2026-11-10_import_ladder/README.md); the larger rate is used (card C-13). Applied to every library. It is
# an estimate for libraries other than Pfam.
SCAN_CPU_MIN_PER_M = 7.96
# Clustering (cards B-36, B-42). Measured, wilkcomp02, 2026-11-10, the clustering run records of the 10 Gbp
# rebuild, 16 threads each (branch build/r10, planning/measurements/2026-11-10_rebuild_10gbp.md, "Per-step numbers,
# phase 1"): `rN.1b` 284.10 CPU-minutes for 2,053,075 proteins and `rN.1` 2,195.91 for 10,620,684. The estimate is
# the power law through those two points (derived):
#     cpu_min(n) = 284.10 × (n / 2,053,075) ** b,   b = ln(2,195.91 / 284.10) / ln(10,620,684 / 2,053,075) = 1.244
# It is exact at both. The third measured point, `rN.1a` (84.75 CPU-minutes for 1,134,550 proteins, phages only, whose
# proteins are shorter), it overestimates by about 60% (derived: 135.8 against 84.75): the safe side for a refusal.
# Above 10.6 M proteins it is extrapolated: no clustering that large has been measured. n is the rows of the release's
# `protein_seq` files, which counts a protein shared by two stages twice (10,634,028 rows for 10,620,684 proteins in
# `rN.1`, measured), so it errs high as well.
CLUSTER_CPU_MIN_LO, CLUSTER_PROTEINS_LO = 284.10, 2_053_075
CLUSTER_CPU_MIN_HI, CLUSTER_PROTEINS_HI = 2195.91, 10_620_684
CLUSTER_EXPONENT = math.log(CLUSTER_CPU_MIN_HI / CLUSTER_CPU_MIN_LO) / math.log(CLUSTER_PROTEINS_HI / CLUSTER_PROTEINS_LO)
# This is the release's estimate for whether clustering runs in parts and for the check of clustering in one process
# (`cluster_refusal`; no exemption there since card B-42). The parts' own cut uses
# the per-call model of `biobox.l2.cluster.call_cpu_min` (card C-22), fitted on the per-call records of the clustering of
# `rN.3` (10.6 M proteins) and `rC.1b` (2.05 M); `tests/r1/test_cluster_call_estimates.py` asserts that the two agree within 20% at both
# measured sizes (the calls and the one-process rest sum to -1.0% and -3.1% of this estimate at 10.6 M and 2.05 M
# proteins, derived). Above 10.6 M both are extrapolations; at 50 M the per-call work (15,772 CPU-minutes, derived)
# is 4.6% above this estimate (15,079).
# The reason the run record of a clustering call run as one process because no number of parts fits it gives (card
# B-45, DEC-29: `biobox.l2.cluster._one_process_doc`). Since card B-42 it is clustering's only reason: the one-process
# clustering step no longer gets DEC-25's `exempt: DEC-25`.
DEC29_EXEMPT = "exempt: DEC-29"
# Scan shards are planned to fill 80% of the CPU-minute budget, not all of it (card C-11): the rate above
# is Pfam's on wilkcomp02, so a shard sized to the whole budget passes it whenever a library is slower
# than Pfam or the host is loaded, and a process over the budget breaks AGENTS §7.2. The margin covers a
# scan up to 25% slower per residue than the rate above (derived: 1 / 0.8 = 1.25; 9.95 CPU-minutes per
# million residues, derived: 7.96 × 1.25) before a shard passes the budget. The same margin decides
# whether a scan is sharded at all.
SHARD_MARGIN = 0.8
# Threads per scan shard when --scan-cpus is not given (or the allowance, if smaller). A default, not a
# measured optimum: the allowance is filled by running more shards at once, not by more threads each.
DEFAULT_SCAN_CPUS = 4
# Seconds between checks of the allowance while a launch is held for it and nothing runs to wait on (card C-12). The
# allowance changes at the edges of the day window and with the load average, so a minute's delay is small beside a
# step or a shard. Not measured; never below a second (no tighter polling).
RECHECK_S = 60.0
# The ncRNA scan (cards E-07, E-08, DEC-28): 304.0 CPU-hours per 10 Gbp, the slower of the two shards of card E-07's
# pod check (*derived* there: shard 1's CPU-seconds / its bases × 1e10 / 3,600, from run records *measured* on
# wilkcomp02, 2026-11-11, planning/measurements/2026-11-11_ncrna_producer.md; shard 0 ran at 206.6, both together at
# 254.5; the design note's per-model derivation gave 247.4). Per-base cost follows the genomes' intergenic records, so
# the slower measured rate is used (card E-08): a shard cut at SHARD_MARGIN of the budget at this rate uses 240 of
# 300 CPU-minutes when it runs at the slower rate, where at 247.4 it would have used 295 (derived: 240 × 304.0 /
# 247.4). Shard size: 131,578,947 genome bases at 300 CPU-minutes (derived: 0.8 × 300 × 1e10 / (304.0 × 60)).
NCRNA_CPU_H_PER_10GBP = 304.0
# Peak RSS of one `cmsearch --cpu 2` on intergenic records, 6.4 GiB at most (*measured*, wilkcomp02, 2026-11-11,
# measurement 0-09c, full library; 4.8 GiB on whole genomes). A shard of C threads is counted at 6.4 GiB per two
# threads (derived; an over-estimate, since card E-01 measured 12.3 GiB at --cpu 8 on whole genomes, not 25.6).
NCRNA_JOB_GIB, NCRNA_JOB_THREADS = 6.4, 2
CUTOFFS = ("gathering", "trusted", "noise", "none")
_READ = 1 << 20  # bytes per read when decompressing an input to count it

Hook = Callable[[str], None]
Runner = Callable[..., int]  # runstep.run_step's signature


@dataclass(frozen=True)
class Library:
    """One domain library to scan: its name and version (the `library=<name>-<version>` partition),
    its HMM file and the bit-score cutoff to apply."""
    name: str
    version: str
    hmm: Path
    cutoff: str


@dataclass(frozen=True)
class Build:
    """The command line of one `build-stage` call."""
    root: Path
    stage: str
    release: str
    genomes: Path
    objects: Path
    stages: tuple[str, ...] | None = None  # None: the release's release_def names them
    libraries: tuple[Library, ...] = ()
    tmp: Path | None = None  # None: <profile tmp>/scratch
    threads: int | None = None  # None: min(max_threads_day, max_threads_night)
    gene_call_version: str = GENE_CALL_VERSION
    scan_jobs: int | None = None  # J, scan shards at once; None: the allowance // C (`scan_slots`)
    scan_cpus: int | None = None  # C, threads per scan shard; None: DEFAULT_SCAN_CPUS, or the allowance if smaller
    shard_bases: int | None = None  # ingest's batch size (`--shard-bases`); None: the store's default (card B-30)
    cluster_parts: int | None = None  # clustering in parts, every per-query call in this many (card B-41); None: auto
    ncrna_cm: Path | None = None  # the Rfam library of the ncRNA producer in the `features` slot (card E-07); None: none
    ncrna_set: str = ncrna.DEFAULT_SET  # its model set: a frozen set's name or a set file


@dataclass(frozen=True)
class Step:
    run_stage: str  # run-step's --stage: the stage, or the release for a release step
    name: str  # run-step's --step
    argv: tuple[str, ...] | None  # the bio-box arguments; None for an empty slot
    check: str | None = None  # "ingest", "scan", "ncrna" or "cluster": the CPU-minute estimate made before it
    module: str = "biobox.cli"  # what `python -m` runs: the bio-box CLI, or biobox.l1.domains for a rescan
    scratch: str | None = None  # the step's own scratch directory (its --tmp, card C-13), or None


def threads_of(prof: Profile, b: Build) -> int:
    return b.threads if b.threads is not None else min(prof.max_threads_day, prof.max_threads_night)


def scan_allowance(prof: Profile, now: datetime, load1: float) -> int:
    """The threads `run-step` allows at local time `now` with this 1-minute load average
    (`runstep.choose_threads`): what the scan shards running at once share."""
    return runstep.choose_threads(prof, now, load1)


def scan_slots(b: Build, allowance: int) -> tuple[int, int]:
    """(J, C): scan shards at once and threads per shard. C is `--scan-cpus`, else `DEFAULT_SCAN_CPUS` or
    the allowance if smaller; J is `--scan-jobs`, else allowance // C, so one call uses the allowance and
    no more. Not checked here (`_refusal` does)."""
    cpus = b.scan_cpus if b.scan_cpus is not None else min(DEFAULT_SCAN_CPUS, allowance)
    jobs = b.scan_jobs if b.scan_jobs is not None else (allowance // cpus if cpus > 0 else 0)
    return jobs, cpus


def scratch_of(prof: Profile, b: Build) -> Path:
    return (Path(b.tmp) if b.tmp is not None else Path(prof.tmp) / "scratch").resolve()


def scan_name(lib: Library) -> str:
    return f"domains-{lib.name}-{lib.version}"


def ncrna_set_arg(b: Build) -> str:
    """`--ncrna-set` as the steps pass it: a set file given by path, made absolute; else the frozen set's name."""
    p = Path(b.ncrna_set)
    return str(p.resolve()) if p.is_file() else b.ncrna_set


def ncrna_name(b: Build) -> str:
    """The ncRNA step's name, `ncrna-<producer>`: a new producer (library, set or region rule) is a new step."""
    return f"ncrna-{ncrna.producer(Path(b.ncrna_cm), ncrna_set_arg(b))}"


def ncrna_cpus(prof: Profile, b: Build) -> int:
    """cmsearch threads of the ncRNA scan: `--scan-cpus`, else `DEFAULT_SCAN_CPUS`, at most `--threads`."""
    return min(threads_of(prof, b), b.scan_cpus if b.scan_cpus is not None else DEFAULT_SCAN_CPUS)


def steps(prof: Profile, b: Build) -> list[Step]:
    """The fixed list of steps for `b` (module docstring)."""
    root, n, tmp = str(Path(b.root).resolve()), str(threads_of(prof, b)), str(scratch_of(prof, b))
    common = ("--root", root, "--gene-call-version", b.gene_call_version)
    stages = () if b.stages is None else ("--stages", ",".join(b.stages))
    prods: tuple[str, ...] = ()
    features = Step(b.stage, "features", None)  # other feature producers: none
    if b.ncrna_cm is not None:
        name = ncrna_name(b)
        prods = ("--producers", f"{b.gene_call_version},{name.removeprefix('ncrna-')}")

    def st(sub: str, *extra: str) -> tuple[str, ...]:
        return ("stage", sub, *common, "--stage", b.stage, *extra)

    def rl(sub: str, *extra: str) -> tuple[str, ...]:
        return ("release", sub, *common, "--release", b.release, *extra)

    def own(run_stage: str, name: str, argv: Callable[[str], tuple[str, ...]], **kw) -> Step:
        """A step given its own scratch directory, `<scratch>/<run stage>.<step>`, as `--tmp` (card C-13)."""
        mine = str(Path(tmp) / f"{run_stage}.{name}")
        return Step(run_stage, name, argv(mine), scratch=mine, **kw)

    s, r = b.stage, b.release
    if b.ncrna_cm is not None:
        features = own(s, name, lambda t: st("ncrna", "--cm", str(Path(b.ncrna_cm).resolve()), "--ncrna-set",
                                                      ncrna_set_arg(b), "--cpus", str(ncrna_cpus(prof, b)), "--tmp", t),
                       check="ncrna")
    return [
        Step(s, "inputs", st("inputs", "--genomes", str(Path(b.genomes).resolve()))),
        own(s, "ingest", lambda t: st("ingest", "--objects", str(Path(b.objects).resolve()), "--tmp", t,
                                      *_shard_bases(b)), check="ingest"),
        own(s, "genecall", lambda t: st("genecall", "--tmp", t, "--workers", n)),
        features,
        own(s, "igr", lambda t: st("igr", "--tmp", t, *prods)),
        own(s, "tapes", lambda t: st("tapes", "--tmp", t, *prods)),
        *(Step(s, scan_name(lib), st("domains", "--library", lib.name, "--library-version", lib.version,
                                     "--hmm", str(Path(lib.hmm).resolve()), "--cutoff", lib.cutoff, "--cpus", n),
               check="scan") for lib in b.libraries),
        own(r, "cluster", lambda t: rl("cluster", *stages, "--tmp", t, "--threads", n), check="cluster"),
        own(r, "tokens", lambda t: rl("tokens", *stages, "--tmp", t, *prods)),
        own(r, "postings", lambda t: rl("postings", *stages, "--tmp", t, *prods)),
        own(r, "dim", lambda t: rl("dim", *stages, "--tmp", t)),
        Step(r, "manifest", rl("manifest")),
        Step(r, "catalog", rl("catalog")),
    ]


def _say(msg: str) -> None:
    runstep.log("build-stage", msg)  # stderr, timestamped, flushed (card C-13)


def _record(root: Path, step: Step) -> Path:
    return root / "runs" / f"stage={step.run_stage}" / f"{step.name}.json"


def fasta_bytes(path: Path) -> int:
    """Bytes of FASTA in `path`, decompressed when its name ends in `.gz` (as `fasta` opens it),
    counted while streaming."""
    path = Path(path)
    if not str(path).endswith(".gz"):
        return path.stat().st_size
    total = 0
    with gzip.open(path, "rb") as fh:
        while chunk := fh.read(_READ):
            total += len(chunk)
    return total


def input_bytes(objects: Path) -> int:
    """Bytes of FASTA that the source objects of `--objects` hold. A relative path is taken from the
    directory of the file, as `bio-box stage ingest` takes it."""
    objects = Path(objects)
    total = 0
    for o in json.loads(objects.read_text()):
        p = Path(o["path"])
        total += fasta_bytes(p if p.is_absolute() else objects.parent / p)
    return total


_SIZES: dict[tuple, list[tuple[Path, int, int]]] = {}  # object_sizes, by the objects' paths, sizes and mtimes


def object_sizes(objects: Path) -> list[tuple[Path, int, int]]:
    """(path, bytes on disk, bytes of FASTA) of each source object of `--objects`, in the file's order (paths as
    `input_bytes` takes them). Decompressing is the cost, so a result is kept for the call's objects while none
    of them changes size or mtime (the refusal and the batch decision read it in one call)."""
    objects = Path(objects)
    paths = []
    for o in json.loads(objects.read_text()):
        p = Path(o["path"])
        paths.append(p if p.is_absolute() else objects.parent / p)
    key = tuple((str(p), p.stat().st_size, p.stat().st_mtime_ns) for p in paths)
    if key not in _SIZES:
        if len(_SIZES) > 8:
            _SIZES.clear()
        _SIZES[key] = [(p, p.stat().st_size, fasta_bytes(p)) for p in paths]
    return _SIZES[key]


def ingest_cpu_min(n_bytes: int) -> float:
    return INGEST_FIXED_CPU_MIN + n_bytes * INGEST_CPU_MIN / INGEST_BP


def _shard_bases(b: Build) -> tuple[str, ...]:
    return () if b.shard_bases is None else ("--shard-bases", str(b.shard_bases))


def stage_residues(root: Path, stage: str, gene_call_version: str) -> int:
    """Residues of the distinct proteins of `stage`, which is what a domain scan reads. They are
    summed from the `length` column of the `protein_seq` files listed in the gene-call step's run
    record (an explicit list, AGENTS §3.9), one row group at a time."""
    rec = json.loads((Path(root) / "runs" / f"stage={stage}" / "genecall.json").read_text())
    prefix = f"l0seq/protein_seq/gene_call_version={gene_call_version}/stage={stage}/"
    total = 0
    for o in rec["outputs"]:
        if not o["path"].startswith(prefix):
            continue
        pf = pq.ParquetFile(Path(root) / o["path"])
        for i in range(pf.num_row_groups):
            total += pc.sum(pf.read_row_group(i, columns=["length"]).column("length")).as_py() or 0
    return total


def scan_cpu_min(residues: int) -> float:
    return residues * SCAN_CPU_MIN_PER_M / 1e6


def cluster_cpu_min(proteins: int) -> float:
    """CPU-minutes of a release's clustering in one process (the derived power law above the constants)."""
    return CLUSTER_CPU_MIN_LO * (proteins / CLUSTER_PROTEINS_LO) ** CLUSTER_EXPONENT


def release_proteins(root: Path, stages: Sequence[str], gene_call_version: str) -> tuple[int | None, str | None]:
    """(rows of the `protein_seq` files of `stages`, None), the files those listed in each stage's gene-call run
    record (an explicit list, AGENTS §3.9) and their counts read from the Parquet footers; or (None, the first stage
    without a gene-call run record)."""
    total = 0
    for s in stages:
        path = Path(root) / "runs" / f"stage={s}" / "genecall.json"
        if not path.exists():
            return None, s
        prefix = f"l0seq/protein_seq/gene_call_version={gene_call_version}/stage={s}/"
        total += sum(pq.read_metadata(Path(root) / o["path"]).num_rows
                     for o in json.loads(path.read_text())["outputs"] if o["path"].startswith(prefix))
    return total, None


def cluster_refusal(prof: Profile, release: str, proteins: int,
                    one_process: Sequence[tuple[str, float]] | None = None) -> str | None:
    """Why clustering `proteins` proteins for `release` may not start (DEC-29, card B-42), or None.

    `one_process` (card B-45) is given when the release clusters in parts: the (call, estimate) of each call its plan
    runs as one process so far (`l2cluster.one_process_calls`, from the published cuts). Each of those is checked as
    before its launch (`one_process_refusal`: refused above the profile's `cluster_cpu_min_max`, and refused when the
    profile has no such key), and nothing else: never the release's estimate in one process, which would refuse a
    large release whose calls all fit the budget in parts (B-45's open question, settled by card B-42). A call cut
    later is checked by its cut (the plan holds the cap, or refuses it without one) and again before its launch.

    Without `one_process` the release clusters in one process, which has no exemption since card B-42 (DEC-29 narrowed
    DEC-25's to the calls no number of parts fits): its estimate over `cpu_min_per_process` is refused whatever
    `cluster_cpu_min_max` says. build-stage never routes such a release there (`_cluster_split` clusters in parts above
    `SHARD_MARGIN` of the budget), so this refusal is a guard, and it names the way to cluster in parts."""
    if one_process is not None:  # in parts: only the calls that run as one process may pass the 300 (DEC-29)
        for name, cpu_min in one_process:
            why = one_process_refusal(prof, release, name, cpu_min)
            if why is not None:
                return why
        return None
    est = cluster_cpu_min(proteins)
    if est <= prof.cpu_min_per_process:
        return None
    return (f"clustering of release {release} over {proteins:,} proteins (the protein_seq rows of its stages) is "
            f"estimated at {est:,.1f} CPU-minutes in one process (a power law through {CLUSTER_CPU_MIN_LO} CPU-min at "
            f"{CLUSTER_PROTEINS_LO:,} and {CLUSTER_CPU_MIN_HI} at {CLUSTER_PROTEINS_HI:,} proteins, measured on "
            f"wilkcomp02, 2026-11-10; extrapolated above {CLUSTER_PROTEINS_HI:,}), over the profile's "
            f"cpu_min_per_process of {prof.cpu_min_per_process} (AGENTS §7.2). Clustering in one process has no "
            f"exemption (card B-42: DEC-29 keeps one only for a call that no number of parts fits, under the profile's "
            f"cluster_cpu_min_max): cluster the release in parts, as build-stage does by itself above "
            f"{SHARD_MARGIN:.0%} of the budget, or with --cluster-parts N")


def one_process_refusal(prof: Profile, release: str, call: str, cpu_min: float) -> str | None:
    """Why `call` of release `release`'s clustering in parts, which its cut runs as one process because no number of
    parts fits it (DEC-29, card B-45), estimated at `cpu_min` CPU-minutes (its work and its own k-mer table), may not
    be launched, or None: only under the profile's `cluster_cpu_min_max`, and not above it."""
    how = (f"{call} of release {release}'s clustering runs as one process (DEC-29: no number of parts fits it, since "
           f"each part would rebuild the target's k-mer table, card B-44), estimated at {cpu_min:,.1f} CPU-minutes")
    cap = prof.cluster_cpu_min_max
    if cap is None:
        return (f"{how}; DEC-29 allows that only under the profile's cluster_cpu_min_max, and this profile has no "
                f"cluster_cpu_min_max (the 300 CPU-minute rule of AGENTS §7.2 holds)")
    if cpu_min > cap:
        return (f"{how}, over the profile's cluster_cpu_min_max of {cap}. Raise cluster_cpu_min_max only with the PI's "
                f"agreement")
    return None


def _refusal(prof: Profile, b: Build, root: Path, plan: list[Step], allowance: int, now: datetime) -> str | None:
    """Why the build may not start (module docstring), or None. `allowance` is `scan_allowance` at `now`."""
    for what, name in (("stage", b.stage), ("release", b.release)):
        if not runstep._NAME.fullmatch(name):
            return f"bad {what} name {name!r}"
    names = [s.name for s in plan]
    for lib in b.libraries:
        if not runstep._NAME.fullmatch(scan_name(lib)) or "/" in lib.name + lib.version:
            return f"bad library name or version {lib.name!r} {lib.version!r}"
        if lib.cutoff not in CUTOFFS:
            return f"library {lib.name}: cutoff {lib.cutoff!r} is not one of {', '.join(CUTOFFS)}"
    if len(set(names)) != len(names):
        return "a library is given twice"
    if b.ncrna_cm is not None:
        ns = next(s for s in plan if s.check == "ncrna")
        if not runstep._NAME.fullmatch(ns.name):
            return f"bad ncRNA step name {ns.name!r}"
        try:
            ncrna.read_set(ncrna_set_arg(b))
        except (FileNotFoundError, ValueError) as e:
            return f"--ncrna-set: {e}"
    if threads_of(prof, b) < 1:
        return f"--threads {threads_of(prof, b)}: must be at least 1"
    if b.cluster_parts is not None and (isinstance(b.cluster_parts, bool) or b.cluster_parts < 1):
        return f"--cluster-parts {b.cluster_parts}: must be at least 1"
    jobs, cpus = scan_slots(b, allowance)
    if jobs < 1 or cpus < 1:
        return f"--scan-jobs {jobs} and --scan-cpus {cpus}: each must be at least 1"
    if jobs * cpus > allowance:
        return (f"--scan-jobs {jobs} × --scan-cpus {cpus} = {jobs * cpus} threads at once, but run-step allows "
                f"{allowance} now ({now:%a %H:%M}: the profile's day cap {prof.max_threads_day} on weekdays from "
                f"{prof.day_start_hour}:00 to {prof.day_end_hour}:00, else the night cap {prof.max_threads_night}, "
                f"and the day cap while the load average is above {runstep.BUSY_LOAD}); lower one of them")
    if threads_of(prof, b) > allowance:
        return (f"--threads {threads_of(prof, b)} is above the {allowance} threads run-step allows now "
                f"({now:%a %H:%M}); give at most {min(prof.max_threads_day, prof.max_threads_night)} (the smaller cap, "
                f"allowed at any hour) or leave it out")
    if scratch_of(prof, b).is_relative_to(root):
        return f"{scratch_of(prof, b)}: the scratch directory may not lie under the data root"
    try:
        if b.stages is None:
            stages = catalog.release_stages(root, b.release)
        else:
            stages = catalog.verify_release_def(root, b.release, list(b.stages))
    except (FileNotFoundError, ValueError, TypeError) as e:
        return str(e)
    if b.stage not in stages:
        return f"release {b.release} lists stages {stages}, not stage {b.stage}"
    ingest = next(s for s in plan if s.check == "ingest")
    if not _record(root, ingest).exists():
        sizes = object_sizes(b.objects)
        n = sum(f for _, _, f in sizes)
        big = max(sizes, key=lambda t: (t[2], str(t[0])), default=None)
        if big is not None and ingest_cpu_min(big[2]) > prof.cpu_min_per_process:
            return (f"ingest of stage {b.stage} reads {n:,} bytes of FASTA, estimated at {ingest_cpu_min(n):,.1f} "
                    f"CPU-minutes in one process ({INGEST_CPU_MIN} CPU-min per {INGEST_BP:,} bp plus "
                    f"{INGEST_FIXED_CPU_MIN} CPU-min to start, measured on wilkcomp02, 2026-11-10). In batch "
                    f"processes (card B-30) each source object is read whole by one process, and {big[0].name} "
                    f"alone holds {big[2]:,} bytes, estimated at {ingest_cpu_min(big[2]):,.1f} CPU-minutes, over "
                    f"the profile's cpu_min_per_process of {prof.cpu_min_per_process} (AGENTS §7.2). Split that "
                    f"object into smaller files")
    cluster = next(s for s in plan if s.check == "cluster")
    if not _record(root, cluster).exists():
        n, _ = release_proteins(root, stages, b.gene_call_version)
        if n is not None:  # every stage is gene-called: the estimate is known now (else at the cluster step)
            return cluster_refusal(prof, b.release, n, _one_process(prof, b, root))
    return None


def _one_process(prof: Profile, b: Build, root: Path) -> list[tuple[str, float]] | None:
    """For `cluster_refusal` (card B-45): when the release clusters in parts, the calls its published cuts run as one
    process (DEC-29); None when it clusters in one process."""
    if not _cluster_split(prof, b, root, quiet=True):
        return None
    return l2cluster.one_process_calls(root, b.release)


def _cluster_gate(prof: Profile, b: Build, root: Path) -> str | None:
    """Before the cluster step's command (card B-36): why it may not start, or None. The release's stages are
    `--stages`, else its `release_def` (`_refusal` has checked both)."""
    stages = list(b.stages) if b.stages is not None else catalog.release_stages(root, b.release)
    n, missing = release_proteins(root, stages, b.gene_call_version)
    if n is None:
        return (f"release {b.release}: stage {missing} has no gene-call run record (runs/stage={missing}/genecall.json), "
                f"so its clustering cannot be estimated before launch (AGENTS §7.2)")
    return cluster_refusal(prof, b.release, n, _one_process(prof, b, root))


def shard_residues(prof: Profile) -> int:
    """Residues per scan shard: what fits in `SHARD_MARGIN` of the profile's `cpu_min_per_process` at the
    measured rate (`SCAN_CPU_MIN_PER_M`), at least one. Nothing else: not the jobs, not the threads."""
    return max(1, int(prof.cpu_min_per_process * SHARD_MARGIN * 1e6 / SCAN_CPU_MIN_PER_M))


def _sharded(prof: Profile, b: Build, root: Path, step: Step) -> bool:
    """Whether the scan `step` runs in shards: its estimate in one process exceeds `SHARD_MARGIN` of the
    budget, or an earlier call planned it in shards (its plan is published)."""
    lib = _library(b, step)
    if (root / domains.scan_dir(b.stage, lib.name, lib.version) / "plan.json").exists():
        return True
    residues = stage_residues(root, b.stage, b.gene_call_version)
    est = scan_cpu_min(residues)
    if est <= SHARD_MARGIN * prof.cpu_min_per_process:
        return False
    _say(f"{b.stage}.{step.name}: the scan of {residues:,} residues is estimated at {est:,.1f} CPU-minutes in one "
         f"process ({SCAN_CPU_MIN_PER_M} CPU-min per million residues, measured for Pfam on wilkcomp02, "
         f"2026-11-10), over {SHARD_MARGIN:.0%} of the profile's cpu_min_per_process of {prof.cpu_min_per_process} "
         f"(AGENTS §7.2); it runs in shards of at most {shard_residues(prof):,} residues, one step each")
    return True


def _library(b: Build, step: Step) -> Library:
    return next(lib for lib in b.libraries if scan_name(lib) == step.name)


def scan_steps(prof: Profile, b: Build, step: Step, n_shards: int | None, cpus: int) -> list[Step]:
    """The steps of the scan `step` in shards: its plan, then (once the plan is published and
    `n_shards` known) one step per shard, each scanning with `cpus` threads, then the merge under the
    scan step's own name, so its run record is the scan's."""
    lib = _library(b, step)
    root, n, tmp = str(Path(b.root).resolve()), str(cpus), str(scratch_of(prof, b))

    def dm(action: str, *extra: str) -> tuple[str, ...]:
        return ("stage", "domains", action, "--root", root, "--gene-call-version", b.gene_call_version,
                "--stage", b.stage, "--library", lib.name, "--library-version", lib.version, *extra)

    hmm = str(Path(lib.hmm).resolve())
    plan = Step(b.stage, f"{step.name}.plan", dm("plan", "--hmm", hmm, "--cutoff", lib.cutoff,
                                                 "--shard-residues", str(shard_residues(prof))))
    if n_shards is None:
        return [plan]
    return [plan,
            *(Step(b.stage, f"{step.name}.shard-{k:05d}", dm("shard", "--hmm", hmm, "--shard", str(k), "--cpus", n,
                                                              "--tmp", tmp)) for k in range(n_shards)),
            Step(b.stage, step.name, dm("merge", "--tmp", tmp))]


def unit_cap() -> int | None:
    """Bytes of the memory cap this process runs under (the `pod_run.sh` unit's), or None (`runstep.detect_cap`)."""
    cap = runstep.detect_cap()
    return None if cap is None else cap["memory_max_bytes"]


def build_stage(prof: Profile, b: Build, *, run: Runner | None = None,
                clock: Callable[[], datetime] = datetime.now, load: Callable[[], float] = runstep._load1,
                mem_cap: Callable[[], int | None] = unit_cap, _hook: Hook | None = None) -> int:
    """Run the steps of `b` in order; return 0, the first step's non-zero exit code, 75 for STOP, or 2
    for a refusal (module docstring). `run` (tests only) stands in for `runstep.run_step`, also for the
    scan shards, which it is then called for from one thread each; without it each shard is its own
    `bio-box run-step` process. `clock` and `load` (tests only) give the time and load average the scan
    allowance is decided at. `_hook` (tests only) is called with `after:<step>` after each step that
    exited 0 or was passed over (scan shards in the order they end), and with `hold:<step>` each time a
    step or shard launch is held for the allowance (card C-12). An exception that `_hook` raises
    is not caught; shards then running are waited for first. The clock and the load are read at the
    start, then before each step without a run record, before each shard launch and at each check
    while a launch is held. `mem_cap` (tests only) gives the memory cap the ncRNA shards share (card E-07)."""
    shard_run = run
    run = runstep.run_step if run is None else run
    root = Path(b.root).resolve()
    try:
        plan = steps(prof, b)
    except (FileNotFoundError, ValueError) as e:  # the ncRNA producer's name needs its library, set and binary
        _say(f"refused before any step: {e}")
        return runstep.EXIT_REFUSED
    if (root / "STOP").exists():
        _say(f"{root / 'STOP'} exists; not starting stage {b.stage}")
        return runstep.EXIT_STOP
    now = clock()
    allowance = scan_allowance(prof, now, load())
    why = _refusal(prof, b, root, plan, allowance, now)
    if why is not None:
        _say(f"refused before any step: {why}")
        return runstep.EXIT_REFUSED
    scratch_of(prof, b).mkdir(parents=True, exist_ok=True)
    n = threads_of(prof, b)
    jobs, cpus = scan_slots(b, allowance)

    def allow() -> int:
        """The threads run-step allows now (read again at every call, card C-12)."""
        return scan_allowance(prof, clock(), load())

    def one(step: Step, **kw: str) -> int | None:
        """Run `step` (STOP checked first, then held while its threads do not fit the allowance); its non-zero exit
        code, else None. `kw` goes on to `run` (`budget_reason` for clustering, card B-36)."""
        if (root / "STOP").exists():
            _say(f"{root / 'STOP'} exists; stopping before {step.run_stage}.{step.name}")
            return runstep.EXIT_STOP
        if not _record(root, step).exists() or _intent(root, step).exists():
            code = _hold_step(root, step, n, allow, _hook)
            if code is not None:
                return code
            if step.scratch is not None:
                _fresh_scratch(step)
        cmd: Sequence[str] = [sys.executable, "-m", step.module, *step.argv]
        code = run(prof, root, step.run_stage, step.name, cmd, threads=n, **kw)
        if code != 0:
            _say(f"{step.run_stage}.{step.name} exited {code}; stopping here; the next call resumes at this step"
                 + ("" if step.scratch is None else f", after removing its scratch {step.scratch}"))
            return code
        if step.scratch is not None:
            shutil.rmtree(step.scratch, ignore_errors=True)  # done: nothing reads a step's scratch after it
        return None

    for step in plan:
        if (root / "STOP").exists():
            _say(f"{root / 'STOP'} exists; stopping before {step.run_stage}.{step.name}")
            return runstep.EXIT_STOP
        if step.argv is not None:
            if step.check == "scan" and not _record(root, step).exists() and _sharded(prof, b, root, step):
                code = _scan_in_shards(prof, b, root, step, one, _hook, jobs, cpus, shard_run, allow)
            elif step.check == "ingest" and not _record(root, step).exists() and _batched(prof, b, root):
                code = _ingest_in_batches(prof, b, root, one, _hook, jobs * cpus, shard_run, allow)
            elif step.check == "cluster":
                why = None if _record(root, step).exists() else _cluster_gate(prof, b, root)
                if why is None and not _record(root, step).exists():
                    why = dim.owner_problem(root, b.release, _release_stages(b, root), tmp=scratch_of(prof, b))
                if why is not None:
                    _say(f"refused before {step.run_stage}.{step.name}: {why}")
                    return runstep.EXIT_REFUSED
                if not _record(root, step).exists() and _cluster_split(prof, b, root):
                    code = _cluster_in_parts(prof, b, root, one, _hook, jobs, cpus, shard_run, allow)
                else:  # in one process, within the budget (card B-42: no exemption)
                    code = one(step)
            elif step.check == "ncrna":
                sharded = not _record(root, step).exists() and _ncrna_sharded(prof, b, root, step)
                nc_cpus, nc_jobs = ncrna_cpus(prof, b), jobs if sharded else 1
                if not _record(root, step).exists() or _intent(root, step).exists():
                    cap = mem_cap()
                    nc_jobs = ncrna_jobs(nc_jobs, nc_cpus, cap)
                    if nc_jobs < 1:
                        _say(f"refused before {step.run_stage}.{step.name}: one cmsearch of {nc_cpus} threads is "
                             f"counted at {ncrna_job_bytes(nc_cpus) / 2**30:.1f} GiB, above the memory cap of "
                             f"{cap / 2**30:.1f} GiB this build runs under; give the unit more memory or lower "
                             "--scan-cpus")
                        return runstep.EXIT_REFUSED
                if sharded:
                    code = _ncrna_in_shards(prof, b, root, step, one, _hook, nc_jobs, nc_cpus, shard_run, allow)
                else:
                    code = one(with_seq_stages(step, root, b.stage))
            elif step.name == "genecall":
                code = one(with_seq_stages(step, root, b.stage))
            else:
                code = one(step)
            if code is not None:
                return code
        if _hook is not None:
            _hook(f"after:{step.name}")
    _say(f"stage {b.stage} and release {b.release}: all {len(plan)} steps done")
    return 0


def with_seq_stages(step: Step, root: Path, stage: str) -> Step:
    """The gene-call step reading the sequences of every stage in the root (card C-20, DEC-26): `--seq-stages`
    is the stage itself, then the others whose `seq_locator` exists when the step starts (`store.stages_on_disk`,
    the stages the stage's ingest compared its sequences against). A contig of the stage's genomes that one of them
    stored first is read from there for training and skipped in the output; a stage that holds none of its contigs
    is never read. The files do not depend on which other stages are listed, since one stage owns each contig."""
    seq = [stage, *(s for s in store.stages_on_disk(root) if s != stage)]
    return dataclasses.replace(step, argv=(*step.argv, "--seq-stages", ",".join(seq)))


def _fresh_scratch(step: Step) -> None:
    """Before an attempt of `step` (card C-13): remove what an earlier attempt left in the step's own scratch
    directory (a killed builder never removes its scratch: 1.2 GB after one stopped clustering step in C-10,
    *measured*), then create it empty."""
    d = Path(step.scratch)
    if d.is_symlink() or (d.exists() and not d.is_dir()):
        d.unlink()
    elif d.exists():
        if any(d.iterdir()):
            _say(f"{step.run_stage}.{step.name}: removing {d}, left by an earlier attempt of this step")
        shutil.rmtree(d)
    d.mkdir(parents=True)


def _hold_step(root: Path, step: Step, need: int, allow: Callable[[], int], hook: Hook | None) -> int | None:
    """Wait until a step of `need` threads, which runs alone, fits the allowance `allow()` (card C-12): None then, or
    75 when a `STOP` file appears meanwhile. While it does not fit, check again every `RECHECK_S` seconds."""
    said = None
    while True:
        a = allow()
        if need <= a:
            return None
        if a != said:
            _say(f"{step.run_stage}.{step.name} needs {need} threads, but run-step allows {a} now; holding it, "
                 f"checking again every {RECHECK_S:g} s")
            said = a
        if hook is not None:
            hook(f"hold:{step.name}")
        threading.Event().wait(RECHECK_S)  # nothing runs to wait on: a timed wait, at least a second
        if (root / "STOP").exists():
            _say(f"{root / 'STOP'} exists; stopping before {step.run_stage}.{step.name}")
            return runstep.EXIT_STOP


def _scan_in_shards(prof: Profile, b: Build, root: Path, step: Step, one: Callable[[Step], int | None],
                    hook: Hook | None, jobs: int, cpus: int, run: Runner | None,
                    allow: Callable[[], int]) -> int | None:
    """The scan `step` as its plan step, then its shard steps, up to `jobs` at once with `cpus` threads each and
    within the allowance `allow()` of the hour (`_run_shards`), then the merge step (`scan_steps`); each is a
    `run-step`, so a call after a kill skips every step already published. Before the merge, a shard whose
    scratch run is missing or fails its record is scanned again (`_rescan_lost`). `run` stands in for the shards'
    run-step processes (tests only). The shards' scratch runs are removed once the merge is published."""
    lib = _library(b, step)
    first = scan_steps(prof, b, step, None, cpus)[0]
    code = one(first)
    if code is not None:
        return code
    if hook is not None:
        hook(f"after:{first.name}")
    shards = len(domains.read_plan(root, b.stage, lib.name, lib.version)["shards"])
    subs = scan_steps(prof, b, step, shards, cpus)
    if run is None:
        run = _process_runner(_profile_file(prof, scratch_of(prof, b)))
    code = _run_shards(prof, root, step, subs[1:-1], jobs, cpus, run, hook, allow)
    if code is not None:
        return code
    code = _rescan_lost(prof, b, root, step, subs[1:-1], jobs, cpus, hook, allow)
    if code is not None:
        return code
    code = one(subs[-1])  # the merge; its hook is the scan step's, after this returns
    if code is not None:
        return code
    run_dir = domains.run_path(scratch_of(prof, b), root, b.stage, lib.name, lib.version, 0).parent
    shutil.rmtree(run_dir, ignore_errors=True)
    return None


def _rescan_lost(prof: Profile, b: Build, root: Path, scan: Step, shards: list[Step], jobs: int, cpus: int,
                 hook: Hook | None, allow: Callable[[], int]) -> int | None:
    """Before the merge (card C-12): every shard has its run record, but scratch may have been cleared or damaged
    since (a reboot clears /tmp), and its run records stop `run-step` from scanning it again. Each shard whose
    scratch run is missing, fails the sha256 its published record holds, or cannot be checked (a record without
    one) is scanned again into scratch (`python -m biobox.l1.domains rescan`, `_rescan_process`), as `_launch` runs
    shards: up to `jobs` at once within the allowance, `cpus` threads each, `STOP` checked before each launch, and
    `after:<shard step>.rescan` to `hook` as each ends. The record is not published again: the new run must
    reproduce it, and a mismatch fails that rescan with an error naming the shard. Return None when every run is
    in place, else the failing code (75 for STOP); the merge does not run then."""
    lib = _library(b, scan)
    scratch = scratch_of(prof, b)
    hmm = str(Path(lib.hmm).resolve())
    todo = []
    for k, s in enumerate(shards):
        why = domains.run_problem(root, b.stage, lib.name, lib.version, k, scratch)
        if why is None:
            continue
        _say(f"{b.stage}.{s.name}: its run record exists, but its run in {scratch} {why}; scanning it again into "
             f"scratch (its published record stays; the new run must reproduce it)")
        if "no run sha256" in why:
            _say(f"{b.stage}.{s.name}: note: its record holds no run sha256, so the new run is accepted when its "
                 f"domain and hit counts equal the record's, and the merge checks its row count")
        todo.append(Step(b.stage, f"{s.name}.rescan",
                         ("rescan", "--root", str(root), "--stage", b.stage, "--library", lib.name,
                          "--library-version", lib.version, "--hmm", hmm, "--shard", str(k), "--cpus", str(cpus),
                          "--tmp", str(scratch)), module="biobox.l1.domains"))
    if not todo:
        return None
    code = _launch(prof, root, todo, jobs, cpus, _rescan_process, hook, allow)
    if code is not None and code != runstep.EXIT_STOP:
        _say(f"{b.stage}.{scan.name}: a shard could not be scanned again to its record (above); the merge does not "
             f"run. A record whose rescan differs is not the scan of these inputs: tell the PI")
    return code


# --- the ncRNA scan in shards (card E-07) --------------------------------------------------------------------


def ncrna_cpu_min(bases: int) -> float:
    """CPU-minutes of the ncRNA scan of `bases` genome bases at `NCRNA_CPU_H_PER_10GBP` (derived)."""
    return bases * NCRNA_CPU_H_PER_10GBP * 60 / 1e10


def ncrna_shard_bases(prof: Profile) -> int:
    """Genome bases per ncRNA scan shard: what fits in `SHARD_MARGIN` of the profile's `cpu_min_per_process` at
    `NCRNA_CPU_H_PER_10GBP`, at least one. Nothing else: not the jobs, not the threads."""
    return max(1, int(SHARD_MARGIN * prof.cpu_min_per_process * 1e10 / (NCRNA_CPU_H_PER_10GBP * 60)))


def ncrna_job_bytes(cpus: int) -> int:
    """Memory one cmsearch of `cpus` threads is counted at: `NCRNA_JOB_GIB` per `NCRNA_JOB_THREADS` threads, rounded
    up to whole pairs (constants above)."""
    return int(math.ceil(cpus / NCRNA_JOB_THREADS) * NCRNA_JOB_GIB * 2**30)


def ncrna_jobs(jobs: int, cpus: int, cap_bytes: int | None) -> int:
    """ncRNA shards at once: `jobs`, and no more than fit `cap_bytes` (the memory cap the build runs under) at
    `ncrna_job_bytes(cpus)` each; `jobs` when there is no cap to see. 0 means not even one fits."""
    return jobs if cap_bytes is None else min(jobs, cap_bytes // ncrna_job_bytes(cpus))


def stage_bases(root: Path, stage: str) -> int:
    """Genome bases the stage's ncRNA scan reads: the `contig_len` of the default variants in the stage's
    `contig_index` files listed in its ingest run record (an explicit list, AGENTS §3.9; a contig another stage
    stored is not in it, DEC-26), one row group at a time."""
    rec = json.loads((Path(root) / "runs" / f"stage={stage}" / "ingest.json").read_text())
    prefix = f"{store.INDEX_DIR}/stage={stage}/"
    total = 0
    for o in rec["outputs"]:
        if not o["path"].startswith(prefix):
            continue
        pf = pq.ParquetFile(Path(root) / o["path"])
        for i in range(pf.num_row_groups):
            t = pf.read_row_group(i, columns=["contig_len", "is_default_variant"])
            total += pc.sum(pc.filter(t.column("contig_len"), t.column("is_default_variant"))).as_py() or 0
    return total


def _ncrna_sharded(prof: Profile, b: Build, root: Path, step: Step) -> bool:
    """Whether the ncRNA scan runs in shards: its estimate in one process exceeds `SHARD_MARGIN` of the budget, or
    an earlier call planned it in shards (its plan is published)."""
    prod = step.name.removeprefix("ncrna-")
    if ncrna.find_plan(root, b.stage, prod) is not None:  # the `ncrna_scan` artifact, or where card E-07 put it
        return True
    bases = stage_bases(root, b.stage)
    est = ncrna_cpu_min(bases)
    if est <= SHARD_MARGIN * prof.cpu_min_per_process:
        return False
    _say(f"{b.stage}.{step.name}: the scan of {bases:,} genome bases is estimated at {est:,.1f} CPU-minutes in one "
         f"process ({NCRNA_CPU_H_PER_10GBP} CPU-hours per 10 Gbp, derived in planning/measurements/"
         "2026-11-11_ncrna_producer.md), "
         f"over {SHARD_MARGIN:.0%} of the profile's cpu_min_per_process of {prof.cpu_min_per_process} (AGENTS §7.2); "
         f"it runs in shards of at most {ncrna_shard_bases(prof):,} bases, one step each")
    return True


def ncrna_scan_steps(prof: Profile, b: Build, step: Step, n_shards: int | None, cpus: int,
                     seq_stages: Sequence[str] = ()) -> list[Step]:
    """The steps of the ncRNA scan `step` in shards: its plan (with its own scratch), then (once the plan is
    published and `n_shards` known) one step per shard, each scanning with `cpus` threads and reading the sequences
    of `seq_stages`, then the merge under the step's own name, so its run record is the scan's."""
    root, tmp = str(Path(b.root).resolve()), str(scratch_of(prof, b))

    def nc(action: str, *extra: str) -> tuple[str, ...]:
        return ("stage", "ncrna", action, "--root", root, "--gene-call-version", b.gene_call_version, "--stage",
                b.stage, "--cm", str(Path(b.ncrna_cm).resolve()), "--ncrna-set", ncrna_set_arg(b), *extra)

    mine = str(Path(tmp) / f"{b.stage}.{step.name}.plan")
    plan = Step(b.stage, f"{step.name}.plan", nc("plan", "--scan-bases", str(ncrna_shard_bases(prof)), "--tmp", mine),
                scratch=mine)
    if n_shards is None:
        return [plan]
    seq = ("--seq-stages", ",".join(seq_stages)) if seq_stages else ()
    return [plan,
            *(Step(b.stage, f"{step.name}.shard-{k:05d}", nc("shard", "--shard", str(k), "--cpus", str(cpus),
                                                             "--tmp", tmp, *seq)) for k in range(n_shards)),
            Step(b.stage, step.name, nc("merge", "--tmp", tmp))]


def _ncrna_in_shards(prof: Profile, b: Build, root: Path, step: Step, one: Callable[[Step], int | None],
                     hook: Hook | None, jobs: int, cpus: int, run: Runner | None,
                     allow: Callable[[], int]) -> int | None:
    """The ncRNA scan `step` as its plan step, its shard steps (up to `jobs` at once, `cpus` threads each, within the
    allowance of the hour, `_run_shards`) and its merge step (`ncrna_scan_steps`), as `_scan_in_shards` runs a domain
    scan; a shard whose scratch runs are lost is scanned again first (`_rescan_lost_ncrna`). The shards' scratch runs
    are removed once the merge is published."""
    prod = step.name.removeprefix("ncrna-")
    first = ncrna_scan_steps(prof, b, step, None, cpus)[0]
    code = one(first)
    if code is not None:
        return code
    if hook is not None:
        hook(f"after:{first.name}")
    shards = len(ncrna.read_plan(root, b.stage, prod)["shards"])
    seq = [b.stage, *(s for s in store.stages_on_disk(root) if s != b.stage)]
    subs = ncrna_scan_steps(prof, b, step, shards, cpus, seq)
    if run is None:
        run = _process_runner(_profile_file(prof, scratch_of(prof, b)))
    code = _run_shards(prof, root, step, subs[1:-1], jobs, cpus, run, hook, allow)
    if code is not None:
        return code
    code = _rescan_lost_ncrna(prof, b, root, step, subs[1:-1], jobs, cpus, seq, hook, allow)
    if code is not None:
        return code
    code = one(subs[-1])  # the merge; its hook is the scan step's, after this returns
    if code is not None:
        return code
    shutil.rmtree(ncrna.run_paths(scratch_of(prof, b), root, b.stage, prod, 0)[0].parent, ignore_errors=True)
    return None


def _rescan_lost_ncrna(prof: Profile, b: Build, root: Path, scan: Step, shards: list[Step], jobs: int, cpus: int,
                       seq: Sequence[str], hook: Hook | None, allow: Callable[[], int]) -> int | None:
    """Before the ncRNA merge: each shard whose scratch runs are missing or fail the sha256 its published record holds
    is scanned again into scratch (`python -m biobox.l1.ncrna rescan`), as `_rescan_lost` does for a domain scan."""
    prod = scan.name.removeprefix("ncrna-")
    scratch = scratch_of(prof, b)
    todo = []
    for k, s in enumerate(shards):
        why = ncrna.run_problem(root, b.stage, prod, k, scratch)
        if why is None:
            continue
        _say(f"{b.stage}.{s.name}: its run record exists, but its runs in {scratch} {why}; scanning it again into "
             f"scratch (its published record stays; the new runs must reproduce it)")
        todo.append(Step(b.stage, f"{s.name}.rescan",
                         ("rescan", "--root", str(root), "--stage", b.stage, "--cm", str(Path(b.ncrna_cm).resolve()),
                          "--ncrna-set", ncrna_set_arg(b), "--shard", str(k), "--cpus", str(cpus), "--tmp",
                          str(scratch), "--seq-stages", ",".join(seq)), module="biobox.l1.ncrna"))
    if not todo:
        return None
    code = _launch(prof, root, todo, jobs, cpus, _rescan_process, hook, allow)
    if code is not None and code != runstep.EXIT_STOP:
        _say(f"{b.stage}.{scan.name}: a shard could not be scanned again to its record (above); the merge does not "
             f"run. A record whose rescan differs is not the scan of these inputs: tell the PI")
    return code


def _rescan_process(prof: Profile, root: Path, stage: str, step: str, cmd: Sequence[str], *, threads: int) -> int:
    """Run a shard rescan `cmd` (a Runner, `_rescan_lost`) as its own process in this process's group, under the
    memory cap and nice a `run-step` would give it (`runstep.memory_cap`, `runstep.build_argv`), with `threads` in
    the thread variables. It is not a run-step: it publishes nothing. A signal gives 128 + its number."""
    try:
        cap = runstep.memory_cap(prof.mem_scope, require=prof.require_mem_cap)
    except runstep.NoCap as e:
        _say(f"not scanning {stage}.{step} again: {e}")
        return runstep.EXIT_REFUSED
    argv = runstep.build_argv(list(cmd), prof.nice_monolithic, prof.mem_scope,
                              inherited=cap if cap["how"] == "inherited" else None)
    env = dict(os.environ, BIOBOX_THREADS=str(threads), **{v: str(threads) for v in runstep.THREAD_VARS})
    code = subprocess.run(argv, env=env, check=False).returncode
    return code if code >= 0 else 128 - code


def _profile_file(prof: Profile, scratch: Path) -> Path:
    """The profile as a YAML file in scratch, named by its sha256, for the shards' `bio-box run-step`
    processes (`profile.load` reads it back to the same values)."""
    data = yaml.safe_dump(dataclasses.asdict(prof), sort_keys=True).encode()
    path = Path(scratch) / f"build-stage-profile-{hashlib.sha256(data).hexdigest()[:16]}.yaml"
    if not path.is_file() or path.read_bytes() != data:
        tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
        tmp.write_bytes(data)
        os.replace(tmp, path)
    return path


def _process_runner(profile_file: Path) -> Runner:
    """A shard runner with `runstep.run_step`'s signature that runs the step as its own `bio-box run-step`
    process in this process's group and waits for it; a process ended by a signal gives 128 + its number."""

    def run(prof: Profile, root: Path, stage: str, step: str, cmd: Sequence[str], *, threads: int) -> int:
        argv = [sys.executable, "-m", "biobox.cli", "run-step", "--profile", str(profile_file), "--root", str(root),
                "--stage", stage, "--step", step, "--threads", str(threads), "--", *cmd]
        code = subprocess.run(argv, check=False).returncode
        return code if code >= 0 else 128 - code

    return run


def _intent(root: Path, step: Step) -> Path:
    return root / "runs" / f"stage={step.run_stage}" / f"{step.name}.publishing.json"


def _run_shards(prof: Profile, root: Path, scan: Step, shards: list[Step], jobs: int, cpus: int, run: Runner,
                hook: Hook | None, allow: Callable[[], int]) -> int | None:
    """Run the shard steps without a run record, up to `jobs` at once (module docstring); return None when
    every shard has its record, else the failed shard's code or 75 for STOP. A shard with a record and no
    intent file is passed over (its hook is called)."""
    todo = []
    for s in shards:
        if _record(root, s).exists() and not _intent(root, s).exists():
            if hook is not None:
                hook(f"after:{s.name}")
        else:
            todo.append(s)
    if not todo:
        return None
    return _launch(prof, root, todo, jobs, cpus, run, hook, allow)


def _launch(prof: Profile, root: Path, todo: list[Step], jobs: int, cpus: int, run: Runner,
            hook: Hook | None, allow: Callable[[], int]) -> int | None:
    """Start the steps of `todo` in order, at most `jobs` running at once, each in a thread that calls `run`
    and puts its result on a queue; wait on the queue (never poll). Before each launch the allowance is read
    again (`allow()`, card C-12): a step is launched only while the steps running and it, at `cpus` threads each,
    fit it; else the launch is held (`hold:<step>` to `hook`), and the queue is waited on with a timeout of
    `RECHECK_S` seconds, after which, or after a step ends, the allowance is read again. A `STOP` file before a
    launch, a step that does not exit 0, or an exception (from a step or `hook`) stops launches; the steps
    running are waited for in every case. An exception is raised again after that."""
    done: queue.SimpleQueue = queue.SimpleQueue()

    def work(k: int, s: Step) -> None:
        try:
            res: int | BaseException = run(prof, root, s.run_stage, s.name,
                                           [sys.executable, "-m", s.module, *s.argv], threads=cpus)
        except BaseException as e:  # handed to the main thread, which raises it once the others have ended
            res = e
        done.put((k, s, res))

    pending = list(enumerate(todo))
    running, failed, error, stopped = 0, {}, None, False
    said = None  # the allowance last named in a hold message
    try:
        while True:
            held = False
            while pending and running < jobs and not (failed or stopped or error is not None):
                if (root / "STOP").exists():
                    _say(f"{root / 'STOP'} exists; launching no more shards, waiting for the {running} running")
                    stopped = True
                    break
                k, s = pending[0]
                a = allow()
                if (running + 1) * cpus > a:
                    held = True
                    if a != said:
                        _say(f"{s.run_stage}.{s.name} needs {cpus} threads beside the {running * cpus} in flight, but "
                             f"run-step allows {a} now; holding it, checking again when a shard ends or every "
                             f"{RECHECK_S:g} s")
                        said = a
                    if hook is not None:
                        hook(f"hold:{s.name}")
                    break
                said = None
                pending.pop(0)
                threading.Thread(target=work, args=(k, s), name=f"build-stage {s.name}").start()
                running += 1
            if not running and not held:
                break
            try:  # until a step's run-step has exited; while a launch is held, at most RECHECK_S seconds
                k, s, res = done.get(timeout=RECHECK_S if held else None)
            except queue.Empty:
                continue
            running -= 1
            if isinstance(res, BaseException):
                error = res if error is None else error
                _say(f"{s.run_stage}.{s.name} raised {res!r}; launching no more shards, waiting for the "
                     f"{running} running")
            elif res != 0:
                failed[k] = res
                _say(f"{s.run_stage}.{s.name} exited {res}; launching no more shards, waiting for the {running} "
                     f"running; the next call runs the shards without a run record"
                     + ("" if s.module == "biobox.cli" else " and scans again each shard whose run is not in scratch"))
            elif hook is not None:
                hook(f"after:{s.name}")
    except BaseException:
        while running:
            done.get()
            running -= 1
        raise
    if error is not None:
        raise error
    if failed:
        return failed[min(failed)]
    return runstep.EXIT_STOP if stopped else None


# --- ingest by batch processes (card B-30) -------------------------------------------------------------------


def _batched(prof: Profile, b: Build, root: Path) -> bool:
    """Whether ingest runs as batch processes: an earlier call published a pass-1 record of the stage, or its
    estimate in one process exceeds `SHARD_MARGIN` of the budget (as a scan is sharded)."""
    if (root / store.plan_dir(b.stage) / "pass1-00000.json").exists():
        return True
    sizes = object_sizes(b.objects)
    n = sum(f for _, _, f in sizes)
    if ingest_cpu_min(n) <= SHARD_MARGIN * prof.cpu_min_per_process:
        return False
    _say(f"{b.stage}.ingest: {n:,} bytes of FASTA are estimated at {ingest_cpu_min(n):,.1f} CPU-minutes in one process "
         f"({INGEST_CPU_MIN} CPU-min per {INGEST_BP:,} bp, measured on wilkcomp02, 2026-11-10), over {SHARD_MARGIN:.0%} "
         f"of the profile's cpu_min_per_process of {prof.cpu_min_per_process} (AGENTS §7.2); it runs as pass-1 shards "
         f"of at most {pass1_bytes(prof, b, root):,} bytes of file, a plan, one step per batch and a merge")
    return True


def pass1_bytes(prof: Profile, b: Build, root: Path) -> int:
    """Bytes of source file per pass-1 shard: the published pass-1 records' (an earlier call cut them), else what
    fits `SHARD_MARGIN` of the budget at the measured rate (`INGEST_CPU_MIN` per `INGEST_BP` bytes of FASTA, after
    the `INGEST_FIXED_CPU_MIN` every process pays to start), in bytes on disk at the stage's ratio of FASTA bytes
    to bytes on disk; at least one. A pure function of the inputs and the profile."""
    if (root / store.plan_dir(b.stage) / "pass1-00000.json").exists():
        return store.read_pass1(root, b.stage, 0)["pass1_bytes"]
    sizes = object_sizes(b.objects)
    on_disk, fasta = sum(d for _, d, _ in sizes), sum(f for _, _, f in sizes)
    ratio = max(1.0, fasta / on_disk) if on_disk else 1.0
    per_shard = SHARD_MARGIN * prof.cpu_min_per_process - INGEST_FIXED_CPU_MIN
    return max(1, int(per_shard * INGEST_BP / INGEST_CPU_MIN / ratio))


def ingest_steps(prof: Profile, b: Build, n_pass1: int, per_shard: int, n_batches: int | None) -> list[Step]:
    """The steps of ingest in batch processes: `ingest.pass1-<k>` for each pass-1 shard, `ingest.plan`, then (once
    the plan is published and `n_batches` known) `ingest.batch-<j>` for each batch and the merge under the step
    name `ingest`, so its run record is the ingest's. Every one is a one-thread process with its own scratch
    directory, `<scratch>/<stage>.<step>`, as `--tmp` (card C-13); the pass-1 runs and the decision table they
    pass on to each other are in `<scratch>` itself (`--runs-tmp`), as a scan's shard runs are."""
    root, tmp, objs = str(Path(b.root).resolve()), str(scratch_of(prof, b)), str(Path(b.objects).resolve())

    def ig(name: str, action: str, *extra: str) -> Step:
        mine = str(Path(tmp) / f"{b.stage}.{name}")
        return Step(b.stage, name, ("stage", "ingest", action, "--root", root, "--gene-call-version",
                                    b.gene_call_version, "--stage", b.stage, "--objects", objs, "--tmp", mine,
                                    "--runs-tmp", tmp, *extra), scratch=mine)

    first = [*(ig(f"ingest.pass1-{k:05d}", "pass1", "--shard", str(k), "--pass1-bytes", str(per_shard),
                  *_shard_bases(b)) for k in range(n_pass1)),
             ig("ingest.plan", "plan", *_shard_bases(b))]
    if n_batches is None:
        return first
    return [*first, *(ig(f"ingest.batch-{j:04d}", "batch", "--batch", str(j)) for j in range(n_batches)),
            ig("ingest", "merge")]


def _own_scratch(run: Runner, todo: list[Step]) -> Runner:
    """`run` around steps launched by `_launch`, each given a fresh own scratch directory before it runs and having
    it removed once it exits 0, as `one` does for a step (card C-13)."""
    mine = {s.name: s for s in todo if s.scratch is not None}

    def wrapped(prof: Profile, root: Path, stage: str, step: str, cmd: Sequence[str], *, threads: int) -> int:
        s = mine.get(step)
        if s is not None:
            _fresh_scratch(s)
        code = run(prof, root, stage, step, cmd, threads=threads)
        if s is not None and code == 0:
            shutil.rmtree(s.scratch, ignore_errors=True)
        return code

    return wrapped


def _ingest_in_batches(prof: Profile, b: Build, root: Path, one: Callable[[Step], int | None], hook: Hook | None,
                       jobs: int, run: Runner | None, allow: Callable[[], int]) -> int | None:
    """Ingest as its pass-1 steps, its plan step, its batch steps and its merge step (`ingest_steps`), each a
    `run-step` of one thread, up to `jobs` at once within the allowance of the hour (`_launch`); a call after a
    kill skips every step already published. A scratch run lost after its record was published is made again
    first (`_rescan_lost_ingest`): before the plan, the pass-1 runs; before the batches, the pass-1 runs and the
    decision table. `run` stands in for the steps' run-step processes (tests only). The scratch runs are removed
    once the merge is published."""
    from biobox.cli import _source_objects

    per_shard = pass1_bytes(prof, b, root)
    # the objects as the steps read them (`--objects`, resolved), so the pass-1 order here is theirs
    n_pass1 = len(store.pass1_shards(_source_objects(str(Path(b.objects).resolve())), per_shard))
    first = ingest_steps(prof, b, n_pass1, per_shard, None)
    plan = first[-1]
    if run is None:
        run = _process_runner(_profile_file(prof, scratch_of(prof, b)))
    jobs = max(1, jobs)
    code = _run_shards(prof, root, plan, first[:-1], jobs, 1, _own_scratch(run, first[:-1]), hook, allow)
    if code is not None:
        return code
    if not _record(root, plan).exists() or _intent(root, plan).exists():
        code = _rescan_lost_ingest(prof, b, root, jobs, hook, allow, plan=False)
        if code is not None:
            return code
    code = one(plan)
    if code is not None:
        return code
    if hook is not None:
        hook(f"after:{plan.name}")
    batches = store.read_ingest_plan(root, b.stage)["batches"]
    for j, x in enumerate(batches):
        if ingest_cpu_min(x["bases_read"]) > prof.cpu_min_per_process:
            _say(f"{b.stage}.ingest.batch-{j:04d}: warning: it reads {x['bases_read']:,} bp, estimated at "
                 f"{ingest_cpu_min(x['bases_read']):,.1f} CPU-minutes, over the budget; a batch is the objects that "
                 f"hold at most shard_bases new bases (DEC-21) and cannot be split")
    subs = ingest_steps(prof, b, n_pass1, per_shard, len(batches))
    code = _rescan_lost_ingest(prof, b, root, jobs, hook, allow, plan=True)
    if code is not None:
        return code
    code = _run_shards(prof, root, plan, subs[n_pass1 + 1:-1], jobs, 1, _own_scratch(run, subs[n_pass1 + 1:-1]), hook,
                       allow)
    if code is not None:
        return code
    code = one(subs[-1])  # the merge; its hook is the ingest step's, after this returns
    if code is not None:
        return code
    shutil.rmtree(store.ingest_runs_dir(scratch_of(prof, b), root, b.stage), ignore_errors=True)
    return None


def _rescan_lost_ingest(prof: Profile, b: Build, root: Path, jobs: int, hook: Hook | None, allow: Callable[[], int],
                        plan: bool) -> int | None:
    """Make again the scratch runs whose record is published but which are missing or fail its sha256
    (`store.ingest_run_problems`; a reboot clears /tmp): each pass-1 run, then (with `plan`) the decision table,
    each by `python -m biobox.l0seq.store rescan` in its own process (`_rescan_process`, as a scan shard's rescan),
    `after:<step>.rescan` to `hook` as each ends. The records stay as published: a rescan must reproduce its
    record's sha256, and one that does not stops the build with an error naming the shard. Return None when
    every run is in place, else the failing code (75 for STOP)."""
    scratch = scratch_of(prof, b)
    lost = store.ingest_run_problems(root, b.stage, scratch)
    base = ("rescan", "--root", str(root), "--stage", b.stage, "--objects", str(Path(b.objects).resolve()),
            "--tmp", str(scratch))
    pass1 = [Step(b.stage, f"ingest.pass1-{k:05d}.rescan", (*base, "--pass1", str(k)), module="biobox.l0seq.store")
             for what, k, _ in lost if what == "pass1"]
    phases = [pass1, [Step(b.stage, "ingest.plan.rescan", (*base, "--plan"), module="biobox.l0seq.store")]
                     if plan and any(what == "plan" for what, _, _ in lost) else []]
    for what, k, why in lost:
        if what == "pass1" or plan:
            _say(f"{b.stage}.ingest.{'plan' if k is None else f'pass1-{k:05d}'}: its record exists, but its run in "
                 f"{scratch} {why}; making it again into scratch (its published record stays; the new run must "
                 f"reproduce it)")
    for todo in phases:
        if not todo:
            continue
        code = _launch(prof, root, todo, jobs, 1, _rescan_process, hook, allow)
        if code is not None:
            if code != runstep.EXIT_STOP:
                _say(f"{b.stage}.ingest: a scratch run could not be made again to its record (above); the steps "
                     f"after it do not run. A record whose rescan differs is not the ingest of these inputs: tell "
                     f"the PI")
            return code
    return None


# --- clustering in parts (card B-41) -------------------------------------------------------------------------


def _release_stages(b: Build, root: Path) -> list[str]:
    return list(b.stages) if b.stages is not None else catalog.release_stages(root, b.release)


def _cluster_split(prof: Profile, b: Build, root: Path, quiet: bool = False) -> bool:
    """Whether the release clusters in parts: an earlier call planned it so (its plan is published), the call forces
    it (`--cluster-parts`), or its estimate in one process exceeds `SHARD_MARGIN` of the budget (as a scan is
    sharded and ingest batched). Every stage of the release is gene-called, so the protein count is known. `quiet`:
    no log line (the gate asks it, card B-45)."""
    if (root / l2cluster.plan_dir(b.release) / f"{l2cluster.PLAN_STEP}.json").exists() or b.cluster_parts is not None:
        return True
    n, _ = release_proteins(root, _release_stages(b, root), b.gene_call_version)
    est = cluster_cpu_min(n)
    if est <= SHARD_MARGIN * prof.cpu_min_per_process:
        return False
    if quiet:
        return True
    _say(f"{b.release}.cluster: clustering {n:,} proteins is estimated at {est:,.1f} CPU-minutes in one process, over "
         f"{SHARD_MARGIN:.0%} of the profile's cpu_min_per_process of {prof.cpu_min_per_process} (AGENTS §7.2); it "
         f"runs in parts (card B-41): a plan, one step per module call, each per-query call cut into parts of at "
         f"most {SHARD_MARGIN * prof.cpu_min_per_process:g} CPU-minutes, then the tables")
    return True


def cluster_steps(prof: Profile, b: Build) -> Callable[..., Step]:
    """A maker of the steps of clustering in parts: `make("plan")`, `make("step", name, threads)`,
    `make("finish")` and `make("redo", name)`, each `bio-box release cluster <action>` with the scratch directory
    (`--tmp`: the work directory is under it, the same for every step) and the plan's cut: `--cluster-parts`, else
    `SHARD_MARGIN` of the budget per part, and the profile's `cluster_cpu_min_max` as the one-process cap when it
    sets one (card B-45, DEC-29)."""
    root, scratch, n = str(Path(b.root).resolve()), str(scratch_of(prof, b)), threads_of(prof, b)
    stages = () if b.stages is None else ("--stages", ",".join(b.stages))
    cut = (("--parts", str(b.cluster_parts)) if b.cluster_parts is not None
           else ("--part-cpu-min", repr(SHARD_MARGIN * prof.cpu_min_per_process)))
    if prof.cluster_cpu_min_max is not None:  # card B-45: a call no number of parts fits runs as one process under it
        cut += ("--one-process-cpu-min-max", repr(float(prof.cluster_cpu_min_max)))

    def cl(action: str, *extra: str) -> tuple[str, ...]:
        return ("release", "cluster", action, "--root", root, "--gene-call-version", b.gene_call_version,
                "--release", b.release, *extra, "--tmp", scratch)

    def make(action: str, name: str | None = None, threads: int | None = None) -> Step:
        if action == "plan":
            return Step(b.release, l2cluster.PLAN_STEP, cl("plan", *stages, *cut))
        if action == "finish":
            return Step(b.release, l2cluster.FINAL_STEP, cl("finish", "--threads", str(n)))
        if action == "redo":
            return Step(b.release, f"{name}.redo", cl("redo", "--name", name, "--threads", str(n)))
        return Step(b.release, name, cl("step", "--name", name, "--threads", str(n if threads is None else threads)))

    return make


def _cluster_in_parts(prof: Profile, b: Build, root: Path, one: Callable[..., int | None], hook: Hook | None,
                      jobs: int, cpus: int, run: Runner | None, allow: Callable[[], int]) -> int | None:
    """The release's clustering in parts (card B-41, `biobox.l2.cluster`): its plan step, then each module call in
    order, a per-query call as its cut step and, when cut into more than one part, its part steps (up to `jobs` at
    once with `cpus` threads each within the allowance of the hour, as scan shards run) and its join, each other
    call as one step; then `cluster` (the tables) under the cluster step's own name, so its run record is the
    clustering's. Every step is a `run-step` of `--threads` threads (the parts of `cpus`), so a call after a kill
    skips each step already published. A call that finds the plan published first makes again, in order, every
    published step whose outputs are not in scratch (`_redo_lost_cluster`). A join's parts are removed once it is
    published; the work directory once the tables are. `run` stands in for the parts' run-step processes (tests)."""
    make = cluster_steps(prof, b)
    scratch = scratch_of(prof, b)
    first = make("plan")
    if _record(root, first).exists() and not _intent(root, first).exists():
        code = _redo_lost_cluster(prof, b, root, hook, allow)
        if code is not None:
            return code
    code = one(first)
    if code is not None:
        return code
    if hook is not None:
        hook(f"after:{first.name}")
    if run is None:
        run = _process_runner(_profile_file(prof, scratch))
    for c in l2cluster.calls():
        parts, kw = 1, {}
        if c.split:
            cut = make("step", f"{c.name}.cut")
            code = one(cut)
            if code is not None:
                return code
            if hook is not None:
                hook(f"after:{cut.name}")
            rec = l2cluster.read_record(root, b.release, cut.name)
            parts = len(rec["parts"])
            if parts > 1:
                subs = [make("step", l2cluster.part_name(c.name, k), cpus) for k in range(parts)]
                code = _run_shards(prof, root, cut, subs, jobs, cpus, run, hook, allow)
                if code is not None:
                    return code
            if "one_process" in rec:  # card B-45: no number of parts fits it; one process under the cap (DEC-29)
                why = one_process_refusal(prof, b.release, c.name, rec["one_process"]["cpu_min"])
                if why is not None and not _record(root, make("step", c.name)).exists():
                    _say(f"refused before {b.release}.{c.name}: {why}")
                    return runstep.EXIT_REFUSED
                kw = {"budget_reason": DEC29_EXEMPT}
        s = make("step", c.name)
        code = one(s, **kw)
        if code is not None:
            return code
        if parts > 1:
            l2cluster.clean_parts(root, b.release, c.name, scratch)
        if hook is not None:
            hook(f"after:{s.name}")
    code = one(make("finish"))  # its hook is the cluster step's, after this returns
    if code is not None:
        return code
    shutil.rmtree(l2cluster.work_dir(scratch, root, b.release), ignore_errors=True)
    return None


def _redo_lost_cluster(prof: Profile, b: Build, root: Path, hook: Hook | None, allow: Callable[[], int]) -> int | None:
    """Before clustering in parts goes on after its plan is published (a reboot clears /tmp): every published step
    whose outputs are not in the work directory as its record holds them (`l2cluster.lost`) is made again, in step
    order, one at a time (`bio-box release cluster redo`, each its own process under the cap and nice of a run-step,
    `_rescan_process`; `STOP` checked before each), `after:<step>.redo` to `hook` as each ends. The records stay as
    published: each must be reproduced, and one that is not stops the build naming the step. Return None when every
    output is in place, else the failing code (75 for STOP)."""
    scratch = scratch_of(prof, b)
    names = l2cluster.lost(root, b.release, scratch)
    if not names:
        return None
    _say(f"{b.release}.cluster: {len(names)} published step(s) have no outputs in scratch as their records hold them, "
         f"from {names[0]} (in {l2cluster.work_dir(scratch, root, b.release)}); making them again in order, each to "
         f"its record")
    make = cluster_steps(prof, b)
    code = _launch(prof, root, [make("redo", x) for x in names], 1, threads_of(prof, b), _rescan_process, hook, allow)
    if code is not None and code != runstep.EXIT_STOP:
        _say(f"{b.release}.cluster: a step could not be made again to its record (above); clustering does not go on. "
             f"A record that is not reproduced is not the clustering of these inputs: tell the PI")
    return code
