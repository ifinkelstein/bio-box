"""`bio-box build-stage`: a stage and its release as a fixed list of `run-step` steps (card C-04).

    bio-box build-stage --profile P --root R --stage S --release REL [--stages S0,S1]
        --genomes genomes.json --objects objects.json [--library NAME VERSION HMM CUTOFF ...]
        [--tmp D] [--threads N] [--gene-call-version V]

No workflow engine: the list `steps` returns and the run records `run-step` writes. Each step is
`runstep.run_step` running `python -m biobox.cli stage ...` or `python -m biobox.cli release ...`
with `--out` left to `$BIOBOX_STAGING` (card C-08), in this order:

    inputs, ingest, genecall, features, igr, tapes, domains-<lib>-<ver> (one per --library),
    cluster, tokens, postings, dim, manifest, catalog

`features` is the slot for feature producers other than the gene caller. It has no command yet
and is passed over. The stage steps are recorded under `runs/stage=<S>/`, the release steps
under `runs/stage=<REL>/`.

The step list and every argument are a pure function of the command line and the profile. Paths
are made absolute. The thread count is `--threads`, else the smaller of the profile's day and
night caps, which `run-step` accepts at any hour; it is given to `run-step` (`BIOBOX_THREADS`), to
gene calling (`--workers`), to the scans (`--cpus`) and to MMseqs2 (`--threads`). Scratch is
`--tmp`, else `<profile tmp>/scratch`. The release steps take the release's stages from its
`release_def`: `--stages` is passed to the builder release steps only when it is given. On first
use it defines the release (the cluster step writes the definition), and later the steps verify it.

Before any step (exit 2, nothing run): a bad name; `--stage` not among the release's stages; no
`--stages` and no `release_def`; `--stages` that disagree with the `release_def`; and ingest whose
estimated CPU-minutes exceed the profile's `cpu_min_per_process` (ingest is one process until card
B-30). The estimate is skipped when ingest's run record exists. Before each step: a `STOP` file in
the root exits 75, as `run-step` does. Before each scan step whose record does not exist: residues
are known only after gene calling, so the scan's estimate is made there, and the build stops before
the scan (exit 2) when it exceeds the budget (one process per library until card B-24). The earlier
steps stay published.

The first step that does not exit 0 ends the build with its exit code. `run-step` gives every step
kill-and-resume (cards C-02, C-03), so the next call skips each step whose run record and outputs
exist and resumes at the step that stopped. A second call after a complete build starts no
command. `build-stage` runs inside the memory-capped unit that `scripts/pod_run.sh` opens on the
host, and every step's `run-step` inherits that cap. It never calls `pod_run.sh` and never polls a
unit.
"""

from __future__ import annotations

import gzip
import json
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq

from biobox import catalog, runstep
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


@dataclass(frozen=True)
class Step:
    run_stage: str  # run-step's --stage: the stage, or the release for a release step
    name: str  # run-step's --step
    argv: tuple[str, ...] | None  # the bio-box arguments; None for an empty slot
    check: str | None = None  # "ingest" or "scan": the CPU-minute estimate made before it


def threads_of(prof: Profile, b: Build) -> int:
    return b.threads if b.threads is not None else min(prof.max_threads_day, prof.max_threads_night)


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
        Step(s, "ingest", st("ingest", "--objects", str(Path(b.objects).resolve()), "--tmp", tmp), check="ingest"),
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


def ingest_cpu_min(n_bytes: int) -> float:
    return n_bytes * INGEST_CPU_MIN / INGEST_BP


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


def _refusal(prof: Profile, b: Build, root: Path, plan: list[Step]) -> str | None:
    """Why the build may not start (module docstring), or None."""
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
        n = input_bytes(b.objects)
        est = ingest_cpu_min(n)
        if est > prof.cpu_min_per_process:
            return (f"ingest of stage {b.stage} reads {n:,} bytes of FASTA, estimated at {est:,.1f} CPU-minutes "
                    f"in one process ({INGEST_CPU_MIN} CPU-min per {INGEST_BP:,} bp, measured on wilkcomp02, "
                    f"2026-11-10), over the profile's cpu_min_per_process of {prof.cpu_min_per_process} "
                    f"(AGENTS §7.2). Split the stage into smaller ones; ingest by batch processes is card B-30")
    return None


def _scan_refusal(prof: Profile, b: Build, root: Path, step: Step) -> str | None:
    residues = stage_residues(root, b.stage, b.gene_call_version)
    est = scan_cpu_min(residues)
    if est <= prof.cpu_min_per_process:
        return None
    return (f"stopped before {b.stage}.{step.name}: the scan of {residues:,} residues is estimated at "
            f"{est:,.1f} CPU-minutes in one process ({SCAN_CPU_MIN_PER_M} CPU-min per million residues, measured "
            f"for Pfam on wilkcomp02, 2026-11-10), over the profile's cpu_min_per_process of "
            f"{prof.cpu_min_per_process} (AGENTS §7.2). The earlier steps are published; a scan in shard "
            f"processes is card B-24")


def build_stage(prof: Profile, b: Build, *, run: Runner | None = None, _hook: Hook | None = None) -> int:
    """Run the steps of `b` in order; return 0, the first step's non-zero exit code, 75 for STOP, or 2
    for a refusal (module docstring). `run` (tests only) stands in for `runstep.run_step`. `_hook` (tests
    only) is called with `after:<step>` after each step that exited 0 or was passed over. An exception
    that `_hook` raises is not caught."""
    run = runstep.run_step if run is None else run
    root = Path(b.root).resolve()
    plan = steps(prof, b)
    if (root / "STOP").exists():
        _say(f"{root / 'STOP'} exists; not starting stage {b.stage}")
        return runstep.EXIT_STOP
    why = _refusal(prof, b, root, plan)
    if why is not None:
        _say(f"refused before any step: {why}")
        return runstep.EXIT_REFUSED
    scratch_of(prof, b).mkdir(parents=True, exist_ok=True)
    n = threads_of(prof, b)
    for step in plan:
        if (root / "STOP").exists():
            _say(f"{root / 'STOP'} exists; stopping before {step.run_stage}.{step.name}")
            return runstep.EXIT_STOP
        if step.argv is not None:
            if step.check == "scan" and not _record(root, step).exists():
                why = _scan_refusal(prof, b, root, step)
                if why is not None:
                    _say(why)
                    return runstep.EXIT_REFUSED
            cmd: Sequence[str] = [sys.executable, "-m", "biobox.cli", *step.argv]
            code = run(prof, root, step.run_stage, step.name, cmd, threads=n)
            if code != 0:
                _say(f"{step.run_stage}.{step.name} exited {code}; stopping here; the next call resumes at "
                     "this step")
                return code
        if _hook is not None:
            _hook(f"after:{step.name}")
    _say(f"stage {b.stage} and release {b.release}: all {len(plan)} steps done")
    return 0
