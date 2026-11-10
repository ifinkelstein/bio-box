"""`bio-box build-stage`: a stage and its release as a fixed list of `run-step` steps (card C-04).

    bio-box build-stage --profile P --root R --stage S --release REL [--stages S0,S1]
        --genomes genomes.json --objects objects.json [--library NAME VERSION HMM CUTOFF ...]
        [--tmp D] [--threads N] [--scan-jobs J] [--scan-cpus C] [--gene-call-version V]

No workflow engine: the list `steps` returns and the run records `run-step` writes. Each step is
`runstep.run_step` running `python -m biobox.cli stage ...` or `python -m biobox.cli release ...`
with `--out` left to `$BIOBOX_STAGING` (card C-08), in this order:

    inputs, ingest, genecall, features, igr, tapes, domains-<lib>-<ver> (one per --library),
    cluster, tokens, postings, dim, manifest, catalog

`features` is the slot for feature producers other than the gene caller. It has no command yet
and is passed over. A domain scan over the CPU budget runs as several steps (below). The stage
steps are recorded under `runs/stage=<S>/`, the release steps under `runs/stage=<REL>/`.

The step list and every argument are a pure function of the command line and the profile. Paths
are made absolute. The thread count is `--threads`, else the smaller of the profile's day and
night caps, which `run-step` accepts at any hour; it is given to `run-step` (`BIOBOX_THREADS`), to
gene calling (`--workers`), to the scans (`--cpus`) and to MMseqs2 (`--threads`). Scratch is
`--tmp`, else `<profile tmp>/scratch`. The release steps take the release's stages from its
`release_def`: `--stages` is passed to the builder release steps only when it is given. On first
use it defines the release (the cluster step writes the definition), and later the steps verify it.

Before any step (exit 2, nothing run): a bad name; `--stage` not among the release's stages; no
`--stages` and no `release_def`; `--stages` that disagree with the `release_def`; a source object
whose own ingest estimate exceeds the profile's `cpu_min_per_process` (one process reads each object
whole, also in batches, card B-30); and `--scan-jobs` × `--scan-cpus`, or `--threads`, above the threads `run-step` allows at
the start of the call (`runstep.choose_threads`: the profile's cap for the hour and weekday, the day
cap while the host is busy). The ingest estimate is skipped when ingest's run record exists. Before
each step: a `STOP` file in the root exits 75, as `run-step` does.

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
(`l0seq/ingest_plan/`) are a build's own, never in a manifest. Lost scratch is made again before the
plan and before the batches (`_rescan_lost_ingest`, as below), and removed after the merge.

Lost scratch (card C-12). A reboot clears `/tmp`; a shard's run record then still exists while its
run in scratch is gone, so `run-step` would never scan it again and the merge would refuse. Before
the merge, every shard's run is checked against its published record (`domains.run_problem`: missing,
another sha256, or a record without one), and exactly those shards are scanned again into scratch
(`python -m biobox.l1.domains rescan`, each its own process under the cap and nice of a run-step,
run like the shards: up to J at once within the allowance, `STOP` before each launch). The records
stay as published: a rescan must reproduce its record's sha256, and one that does not stops the
build before the merge with an error naming the shard.

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

from biobox import catalog, runstep
from biobox.l0seq import store
from biobox.l1 import domains
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.profile import Profile

# Measured, wilkcomp02, 2026-11-10: one ingest process took 42.9 CPU-minutes for 9,596,180,328 bp of
# FASTA (stage P2 of the import ladder, planning/measurements/2026-11-10_import_ladder/README.md,
# "Step timings" and finding 4). The estimate applies the rate to bytes of decompressed FASTA, which
# are at least the bp they hold (headers and newlines add to them). The estimate therefore errs high.
INGEST_CPU_MIN = 42.9
INGEST_BP = 9_596_180_328
# Measured, wilkcomp02, 2026-11-10: the Pfam scan of the 918,525 proteins new in stage P1 (294.5 M
# residues) took 2,149.9 CPU-minutes, 7.30 per million residues (same README). Applied to every
# library. It is an estimate for libraries other than Pfam.
SCAN_CPU_MIN_PER_M = 7.30
# Scan shards are planned to fill 80% of the CPU-minute budget, not all of it (card C-11): the rate above
# is Pfam's on wilkcomp02, so a shard sized to the whole budget passes it whenever a library is slower
# than Pfam or the host is loaded, and a process over the budget breaks AGENTS §7.2. The same margin
# decides whether a scan is sharded at all.
SHARD_MARGIN = 0.8
# Threads per scan shard when --scan-cpus is not given (or the allowance, if smaller). A default, not a
# measured optimum: the allowance is filled by running more shards at once, not by more threads each.
DEFAULT_SCAN_CPUS = 4
# Seconds between checks of the allowance while a launch is held for it and nothing runs to wait on (card C-12). The
# allowance changes at the edges of the day window and with the load average, so a minute's delay is small beside a
# step or a shard. Not measured; never below a second (no tighter polling).
RECHECK_S = 60.0
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


@dataclass(frozen=True)
class Step:
    run_stage: str  # run-step's --stage: the stage, or the release for a release step
    name: str  # run-step's --step
    argv: tuple[str, ...] | None  # the bio-box arguments; None for an empty slot
    check: str | None = None  # "ingest" or "scan": the CPU-minute estimate made before it
    module: str = "biobox.cli"  # what `python -m` runs: the bio-box CLI, or biobox.l1.domains for a rescan


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


def steps(prof: Profile, b: Build) -> list[Step]:
    """The fixed list of steps for `b` (module docstring)."""
    root, n, tmp = str(Path(b.root).resolve()), str(threads_of(prof, b)), str(scratch_of(prof, b))
    common = ("--root", root, "--gene-call-version", b.gene_call_version)
    stages = () if b.stages is None else ("--stages", ",".join(b.stages))

    def st(sub: str, *extra: str) -> tuple[str, ...]:
        return ("stage", sub, *common, "--stage", b.stage, *extra)

    def rl(sub: str, *extra: str) -> tuple[str, ...]:
        return ("release", sub, *common, "--release", b.release, *extra)

    s, r = b.stage, b.release
    return [
        Step(s, "inputs", st("inputs", "--genomes", str(Path(b.genomes).resolve()))),
        Step(s, "ingest", st("ingest", "--objects", str(Path(b.objects).resolve()), "--tmp", tmp, *_shard_bases(b)),
             check="ingest"),
        Step(s, "genecall", st("genecall", "--tmp", tmp, "--workers", n)),
        Step(s, "features", None),  # other feature producers: none yet
        Step(s, "igr", st("igr", "--tmp", tmp)),
        Step(s, "tapes", st("tapes", "--tmp", tmp)),
        *(Step(s, scan_name(lib), st("domains", "--library", lib.name, "--library-version", lib.version,
                                     "--hmm", str(Path(lib.hmm).resolve()), "--cutoff", lib.cutoff, "--cpus", n),
               check="scan") for lib in b.libraries),
        Step(r, "cluster", rl("cluster", *stages, "--tmp", tmp, "--threads", n)),
        Step(r, "tokens", rl("tokens", *stages, "--tmp", tmp)),
        Step(r, "postings", rl("postings", *stages, "--tmp", tmp)),
        Step(r, "dim", rl("dim", *stages, "--tmp", tmp)),
        Step(r, "manifest", rl("manifest")),
        Step(r, "catalog", rl("catalog")),
    ]


def _say(msg: str) -> None:
    print(f"build-stage: {msg}", file=sys.stderr)


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
    return n_bytes * INGEST_CPU_MIN / INGEST_BP


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
    if threads_of(prof, b) < 1:
        return f"--threads {threads_of(prof, b)}: must be at least 1"
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
                    f"CPU-minutes in one process ({INGEST_CPU_MIN} CPU-min per {INGEST_BP:,} bp, measured on "
                    f"wilkcomp02, 2026-11-10). In batch processes (card B-30) each source object is read whole by "
                    f"one process, and {big[0].name} alone holds {big[2]:,} bytes, estimated at "
                    f"{ingest_cpu_min(big[2]):,.1f} CPU-minutes, over the profile's cpu_min_per_process of "
                    f"{prof.cpu_min_per_process} (AGENTS §7.2). Split that object into smaller files")
    return None


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


def build_stage(prof: Profile, b: Build, *, run: Runner | None = None,
                clock: Callable[[], datetime] = datetime.now, load: Callable[[], float] = runstep._load1,
                _hook: Hook | None = None) -> int:
    """Run the steps of `b` in order; return 0, the first step's non-zero exit code, 75 for STOP, or 2
    for a refusal (module docstring). `run` (tests only) stands in for `runstep.run_step`, also for the
    scan shards, which it is then called for from one thread each; without it each shard is its own
    `bio-box run-step` process. `clock` and `load` (tests only) give the time and load average the scan
    allowance is decided at. `_hook` (tests only) is called with `after:<step>` after each step that
    exited 0 or was passed over (scan shards in the order they end), and with `hold:<step>` each time a
    step or shard launch is held for the allowance (card C-12). An exception that `_hook` raises
    is not caught; shards then running are waited for first. The clock and the load are read at the
    start, then before each step without a run record, before each shard launch and at each check
    while a launch is held."""
    shard_run = run
    run = runstep.run_step if run is None else run
    root = Path(b.root).resolve()
    plan = steps(prof, b)
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

    def one(step: Step) -> int | None:
        """Run `step` (STOP checked first, then held while its threads do not fit the allowance); its non-zero exit
        code, else None."""
        if (root / "STOP").exists():
            _say(f"{root / 'STOP'} exists; stopping before {step.run_stage}.{step.name}")
            return runstep.EXIT_STOP
        if not _record(root, step).exists() or _intent(root, step).exists():
            code = _hold_step(root, step, n, allow, _hook)
            if code is not None:
                return code
        cmd: Sequence[str] = [sys.executable, "-m", step.module, *step.argv]
        code = run(prof, root, step.run_stage, step.name, cmd, threads=n)
        if code != 0:
            _say(f"{step.run_stage}.{step.name} exited {code}; stopping here; the next call resumes at this step")
            return code
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
            else:
                code = one(step)
            if code is not None:
                return code
        if _hook is not None:
            _hook(f"after:{step.name}")
    _say(f"stage {b.stage} and release {b.release}: all {len(plan)} steps done")
    return 0


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
    fits `SHARD_MARGIN` of the budget at the measured rate (`INGEST_CPU_MIN` per `INGEST_BP` bytes of FASTA), in
    bytes on disk at the stage's ratio of FASTA bytes to bytes on disk; at least one. A pure function of the
    inputs and the profile."""
    if (root / store.plan_dir(b.stage) / "pass1-00000.json").exists():
        return store.read_pass1(root, b.stage, 0)["pass1_bytes"]
    sizes = object_sizes(b.objects)
    on_disk, fasta = sum(d for _, d, _ in sizes), sum(f for _, _, f in sizes)
    ratio = max(1.0, fasta / on_disk) if on_disk else 1.0
    return max(1, int(SHARD_MARGIN * prof.cpu_min_per_process * INGEST_BP / INGEST_CPU_MIN / ratio))


def ingest_steps(prof: Profile, b: Build, n_pass1: int, per_shard: int, n_batches: int | None) -> list[Step]:
    """The steps of ingest in batch processes: `ingest.pass1-<k>` for each pass-1 shard, `ingest.plan`, then (once
    the plan is published and `n_batches` known) `ingest.batch-<j>` for each batch and the merge under the step
    name `ingest`, so its run record is the ingest's. Every one is a one-thread process."""
    root, tmp, objs = str(Path(b.root).resolve()), str(scratch_of(prof, b)), str(Path(b.objects).resolve())

    def ig(action: str, *extra: str) -> tuple[str, ...]:
        return ("stage", "ingest", action, "--root", root, "--gene-call-version", b.gene_call_version,
                "--stage", b.stage, "--objects", objs, "--tmp", tmp, *extra)

    first = [*(Step(b.stage, f"ingest.pass1-{k:05d}", ig("pass1", "--shard", str(k), "--pass1-bytes", str(per_shard),
                                                          *_shard_bases(b))) for k in range(n_pass1)),
             Step(b.stage, "ingest.plan", ig("plan", *_shard_bases(b)))]
    if n_batches is None:
        return first
    return [*first, *(Step(b.stage, f"ingest.batch-{j:04d}", ig("batch", "--batch", str(j))) for j in range(n_batches)),
            Step(b.stage, "ingest", ig("merge"))]


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
    code = _run_shards(prof, root, plan, first[:-1], jobs, 1, run, hook, allow)
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
    code = _run_shards(prof, root, plan, subs[n_pass1 + 1:-1], jobs, 1, run, hook, allow)
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
