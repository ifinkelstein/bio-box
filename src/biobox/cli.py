"""`bio-box` command line: search pack, search runs, neighborhoods (SEARCH_PLAN §3.3).

    bio-box pack build --root R --release rS0.1 [--model pfam-38.2=/path/Pfam-A.hmm ...] [--pack-dir P]
    bio-box pack verify --root R --release rS0.1 [--pack-dir P]
    bio-box search run recipe.yaml --root R --release rS0.1 [--wait] [--out DIR | --runs DIR] [--pack-dir P]
    bio-box search status|cancel <run> --root R [--runs DIR]
    bio-box search hits <run> --root R [--runs DIR] [--round N] [--members] [--json] [--top_n 50]
    bio-box search hits <run> --root R [--runs DIR] --format tsv|faa [--round N] [--library L] [--out FILE]
    bio-box search hmm --root R --release R1 --hmm FILE --bits 30 [--neighbors 5]
                --positive ACC:START-END|ACC:START..END|TEXT --out DIR [--model NAME]
                [--profile P | --pack-dir P --runs DIR --work DIR]
    bio-box find --root R --release R1 (--accession ACC[:START-END|:START..END] | --protein UID) [--tape-version T]
    bio-box releases --root R [--profile P | --pack-dir P] [--json]
    bio-box seq --root R --release R1 (--protein UID [UID ...] | --proteins-file F
                | --contig ACC[:START-END|:START..END] [--strand -1]) [--out FILE]
    bio-box sql --root R --release R1 (FILE | -e QUERY) [--runs DIR ...] [--profile P] [--format tsv|json]
                [--out FILE]
    bio-box nbr --root R --anchors <run|uid,uid,...> --window 10g --strand any --edge truncate [--json]
                [--tape-version T] [--runs DIR]
    bio-box nbr [summarize] --anchors-file FILE ...   (protein uids one per line, in place of the anchors argument)
    bio-box nbr <run|uid,uid,...> --root R --window 5g --format tsv [--release R] [--out DIR] [--runs DIR]
    bio-box nbr summarize <run|uid,uid,...> --root R --window 5g --by pfam,nr50 [--top 50] [--release R]
                [--json | --format tsv] [--out FILE] [--runs DIR]
    bio-box nbr flanks <run|uid,uid,...> --root R --out DIR [--window 5g] [--level nr50|nr90|pfam]
                [--group-by anchor-nr90|anchor-nr50|none] [--max-loci 2000] [--cut 0.5] [--distance j|a|combined]
                [--export gff] [--figure [--figure-rows 60]] [--release R] [--runs DIR]
                (card W-25, biobox.search.flanks; --figure: card W-27)
    bio-box export proteins --root R --release R [--out FILE]
    bio-box run-step --profile P --root R --stage S --step NAME [--nice monolithic|fanout] [--threads N] -- CMD ...
    bio-box release define --root R --release R1 --stages S0,S1 [--out DIR]
    bio-box release cluster|tokens|postings|dim|manifest|catalog --root R --release R1 [--stages S0,S1]
                [--tmp D] [--threads N] [--producers P,Q] [--libraries L,M] [--tape-version T] [--out DIR]
    bio-box release cluster plan --root R --release R1 [--stages S0,S1] --tmp D (--parts N | --part-cpu-min M)
        [--one-process-cpu-min-max C]
    bio-box release cluster step|redo --root R --release R1 --name STEP --tmp D --threads N
    bio-box release cluster finish --root R --release R1 --tmp D --threads N
    bio-box release verify --root R --release R2 [--previous R1]
    bio-box stage inputs --root R --stage S --genomes genomes.json [--out DIR]
    bio-box stage ingest --root R --stage S --objects objects.json [--shard-bases N] [--line-width N] [--tmp D]
    bio-box stage ingest pass1 --root R --stage S --objects F --shard K --pass1-bytes N [--shard-bases N] --tmp D
    bio-box stage ingest plan --root R --stage S --objects F [--shard-bases N] [--line-width N] --tmp D
    bio-box stage ingest batch --root R --stage S --objects F --batch J --tmp D
    bio-box stage ingest merge --root R --stage S --objects F --tmp D
    bio-box stage genecall|igr|tapes --root R --stage S [--tmp D] [--producers P,Q] [--workers N] [--out DIR]
    bio-box stage domains --root R --stage S --library pfam --library-version 38.2 --hmm F --cutoff gathering
    bio-box stage domains plan --root R --stage S --library L --library-version V --hmm F --cutoff C --shard-residues N
    bio-box stage domains shard --root R --stage S --library L --library-version V --hmm F --shard K --cpus N --tmp D
    bio-box stage domains merge --root R --stage S --library L --library-version V --tmp D
    bio-box stage ncrna --root R --stage S --cm Rfam.cm.gz [--ncrna-set prok-2] --cpus N [--tmp D] [--seq-stages S0,S1]
                [--out DIR]
    bio-box stage ncrna plan --root R --stage S --cm F [--ncrna-set SET] --scan-bases N [--tmp D]
    bio-box stage ncrna shard --root R --stage S --cm F [--ncrna-set SET] --shard K --cpus N --tmp D [--seq-stages S0,S1]
    bio-box stage ncrna merge --root R --stage S --cm F [--ncrna-set SET] --tmp D
    bio-box build-stage --profile P --root R --stage S --release R1 [--stages S0,S1] --genomes genomes.json
                --objects objects.json [--library NAME VERSION HMM CUTOFF ...] [--tmp D] [--threads N]
                [--scan-jobs J] [--scan-cpus C] [--cluster-parts N] [--ncrna-cm Rfam.cm.gz [--ncrna-set prok-2]]

`build-stage` (card C-04, `biobox.stage`) runs every `stage` step and then every `release` step of
one stage and its release, each as a `run-step`; a second call skips what is done.

`release cluster plan|step|finish|redo` (card B-41, `biobox.l2.cluster`) is a release's clustering in
parts, each step its own process: `build-stage` runs them when clustering in one process is estimated
over 80% of the CPU-minute budget. `redo` writes nothing into the root and is not a run-step.

A release names its stages (card B-25): `--stages` defines it once (its `release_def`); after that
every release step reads the stages from `release_def`, and `--stages`, if given, must repeat them.

Every builder step reads its inputs from `--root` and writes its outputs under `--out` (card C-08),
which defaults to `$BIOBOX_STAGING` when that is set and to the root otherwise; so the same command
runs directly or as `bio-box run-step ... -- python -m biobox.cli stage genecall --root R --stage S`.
A step whose partitions already hold files in the root or in `--out` is refused (exit 2) before any work.

The search path writes nothing under the data root (card W-03). The search pack is built into a pack
directory (`--pack-dir`, `biobox.search.pack`); a run into `--out DIR`, else `run=<run_id>` in a runs
directory (`--runs`). A `<run>` argument is a run directory, or a run id in the runs directory. Defaults,
with a host profile (`--profile`, else $BIOBOX_PROFILE): `<tmp>/packs` and `<tmp>/search/<user>/runs`;
without one, `biobox.search.runner.default_pack_base()` and `default_runs_dir()`. A search whose pack is
missing exits 2 with the command that builds it. `search hits` and `nbr` on a closed run whose pack the pack
directory lacks exit 2 naming the pack (card W-24): a run's genes are named in its pack's tape version, which
only the pack's manifest records. `nbr` reads a run with its own release's files; `--release` naming another
release is refused (card W-24), and the run's proteins go to another release as protein uids.

The profile's `<tmp>/packs` is the host's one pack directory (card W-22): `pack`, `search`, `nbr` and
`releases` find a release's pack there with only `--profile`; `--pack-dir` names another. Only a `pack build`
that makes a pack writes into it, so it may be read-only for everyone but its builder: `pack verify` and a
no-op `pack build` sort under the profile's `<tmp>/search/<user>/work` (without a profile, the system
temporary directory), and a build that cannot write the pack directory is refused (exit 2).

`search run --wait` runs the search in this process and writes the same pending record as a submitted
run, so `search status` and `search cancel` work on it from another shell (card W-13). Without `--wait`
the run is submitted: a detached runner, and the command returns at once. Inside a `scripts/pod_run.sh`
unit (BIOBOX_POD_UNIT, which the script sets) that is refused (exit 2): the unit is a transient service
started with `systemd-run --collect`, which ends when its main process exits and stops every process
left in its cgroup, the detached runner included; a unit stays alive only as long as its main process
runs, and with `--wait` the main process is the runner.

`search hmm` (card W-07) is the discovery workflow in one command: it writes a recipe for a user's HMM file
(`hmm_file` seed, `hmmsearch`, one round at `--bits`, window `--neighbors` genes), runs it in this process
as `search run --wait` does (so it may run inside a `pod_run.sh` unit), then writes into `--out DIR` what the
separate commands write: `recipe.yaml` (the recipe it ran), `hits.tsv` and `hits.faa` (`search hits --format
tsv|faa`), `neighbors.tsv` (`nbr summarize --format tsv`, by `pfam,nr50`; by `nr50,nr90` when the release has
no domain library) and `loci/` (`nbr --format tsv`), each by calling that command. DIR must be new or empty
and outside the data root; it is filled under a temporary name beside it and renamed into place. The positive
control (`--positive`) is required. It names one gene as `ACCESSION:START-END` (the accession, start and end
`search hits --format tsv` prints: 0-based, end exclusive), as `ACCESSION:START..END` (GenBank) or is a substring of a contig header; it is resolved
through the release's contig labels and the slots of the pack's tape version, and anything but exactly one
protein is refused (exit 2) with a list of what matched. A missing pack is refused with the `pack build`
command, as `search run` refuses it; it is never built here.

Coordinates (card W-23, `parse_span`): wherever a gene or region is named by accession, `ACC:START-END` is
bio-box's own form (0-based, end exclusive, as every table prints it) and `ACC:START..END` is GenBank's
(1-based, both ends included): `NC_001604:3171..5822` is `NC_001604:3170-5822`. Anything else after the last
colon made of digits, dots and dashes is refused (exit 2) with both spellings. `find` prints both forms.

`releases`, `seq` and `sql` (card W-23) read a release without exporting it. `releases` lists every release
the root's manifests and release definitions name, from the manifest's file list and Parquet footers (genomes
take one column of `contig_dim`: no footer holds a distinct count). `seq` prints proteins by uid, or a contig
slice by accession or contig uid, read by key (`export.protein_seqs`, `contig_slice`). `sql` runs SELECT
statements on the release's catalog (`sql_connection`: read-only; see there for the protections it relies on)
and prints to stderr the row groups each statement read.

Every command that writes an output file given by `--out` (and `search hmm`, for each file in its DIR) prints
to stderr one line per file: `wrote PATH: N rows|sequences|records` (`N tables, M rows` for a directory of
locus tables), card W-23.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import duckdb
import pyarrow.parquet as pq
import yaml

from biobox import catalog
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.search import export, nbr, pack, runner
from biobox.search.recipe import RecipeError


def _print(lines) -> None:
    for line in lines:
        sys.stdout.write(line + "\n")


_WRITTEN: dict[str, str] = {}  # card W-23: output path -> what it holds ("12 rows"), for every file written here
_QUIET = False  # set by `search hmm` while the commands it calls write into its temporary directory


def _wrote(path, text: str) -> None:
    """One stderr line for an output file written (card W-23): `wrote PATH: N rows|sequences|records`."""
    _WRITTEN[str(path)] = text
    if not _QUIET:
        print(f"wrote {path}: {text}", file=sys.stderr)


class _Tally:
    """The lines of an output, counted as they pass (card W-23): `sequences` counts FASTA records, `records` NDJSON
    lines, `rows` the rows of a TSV of blank-line-separated tables, each a column line and its rows (`#` lines are
    not rows)."""

    def __init__(self, lines, unit: str):
        self.lines, self.unit, self.n = lines, unit, 0

    def __iter__(self):
        fresh = True
        for line in self.lines:
            if self.unit == "sequences":
                self.n += line.startswith(">")
            elif self.unit == "records":
                self.n += 1
            elif not line:
                fresh = True
            elif not line.startswith("#"):
                self.n += not fresh
                fresh = False
            yield line


def _emit(lines, out: str | None, unit: str = "rows") -> None:
    """To stdout, or to the file `out` (written under a temporary name and renamed into place), then one line on
    stderr with what it holds, counted in `unit` (`_Tally`)."""
    if out is None:
        _print(lines)
    else:
        t = _Tally(lines, unit)
        export.write_lines(t, Path(out))
        _wrote(out, f"{t.n} {unit}")


def _profile_tmp(a) -> Path | None:
    """The host profile's `tmp` (`--profile`, else $BIOBOX_PROFILE), or None without a profile."""
    path = getattr(a, "profile", None) or os.environ.get("BIOBOX_PROFILE")
    if not path:
        return None
    from biobox import profile

    return Path(profile.load(path).tmp)


def _packs(a) -> Path:
    """The pack directory (card W-03): `--pack-dir`, else the default for the host profile."""
    return Path(a.pack_dir) if getattr(a, "pack_dir", None) else runner.default_pack_base(_profile_tmp(a))


def _runs(a) -> Path:
    """The runs directory (card W-03): `--runs`, else the per-user default for the host profile."""
    return Path(a.runs) if getattr(a, "runs", None) else runner.default_runs_dir(_profile_tmp(a))


def _run_path(a, ref: str) -> Path:
    """A run named on the command line: a run id (32 hex digits) in the runs directory (`_runs`), unless a
    directory of that name exists here; anything else is the run's directory."""
    if len(ref) == 32 and all(c in "0123456789abcdef" for c in ref) and not Path(ref).is_dir():
        return runner.run_dir(_runs(a), ref)
    return Path(ref)


def _run_release(d: Path) -> str | None:
    """The release a recorded run ran against; None when the run is not on disk."""
    path = Path(d) / "run.parquet"
    if not path.exists():
        return None
    return pq.read_table(path, columns=["ran_against"]).column("ran_against")[0].as_py()


class RunPackMissing(ValueError):
    """A closed run's pack is not in the pack directory, so the tape version that names its genes is unknown."""


def _run_tape_version(d: Path, packs: Path) -> str:
    """The tape version a closed run read: its pack's (card B-16), the pack found by the run's pack_id in the pack
    directory `packs`. RunPackMissing, naming the run's release, its pack id and where it was looked for, when that
    pack is not there (card W-24): the gene-call producer's tape version, the fallback before, names other slots once
    a release has a second producer, and the answer was wrong with exit 0."""
    row = pq.read_table(Path(d) / "run.parquet", columns=["run_id", "ran_against", "pack_id"]).to_pylist()[0]
    pdir = pack.pack_dir(packs, row["ran_against"], row["pack_id"])
    if not (pdir / pack.MANIFEST).is_file():
        raise RunPackMissing(
            f"run {row['run_id'].hex()} searched pack {row['pack_id'].hex()} of release {row['ran_against']}, which is "
            f"not in the pack directory {packs} (no {pdir / pack.MANIFEST}); the run's genes are named in that pack's "
            "tape version, which only its manifest records. Give the pack directory that holds it with --pack-dir")
    return pack.manifest_tape_version(pack.read_manifest(pdir)[0])


def _verify_scratch(a) -> Path | None:
    """Where `pack verify` and a no-op `pack build` sort (card W-22): the profile's per-user work directory
    (`<tmp>/search/<user>/work`, where a search's own verify sorts), else None (the system temporary directory);
    never the pack directory, which may be read-only for all but its builder."""
    tmp = _profile_tmp(a)
    return runner.default_work_base(tmp) if tmp is not None else None


def _pack(a) -> int:
    """`pack build|verify` (card W-03): the release's pack in the pack directory, never under the root. Only a
    build that makes a pack writes into the pack directory (card W-22)."""
    packs = _packs(a)
    scratch = _verify_scratch(a)
    if a.action == "build":
        models = dict(m.split("=", 1) for m in a.model)
        try:
            man = pack.build(Path(a.root), a.release, a.gene_call_version, {k: Path(v) for k, v in models.items()},
                             packs=packs, scratch=scratch)
        except (ValueError, FileExistsError) as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        except PermissionError as e:  # card W-22: a pack directory may be read-only for all but its builder
            print(f"refused: {e}. A pack directory may be read-only for all but its builder: ask the owner of "
                  f"{packs} to build the pack, or build it into a pack directory of your own with --pack-dir",
                  file=sys.stderr)
            return 2
        recs = (export.record("pack_target", {"release_id": a.release, "pack_id": man["pack_id"], "target": t,
                                              "shards": sum(1 for s in man["shards"] if s["target"] == t), **info})
                for t, info in sorted(man["targets"].items()))
        _print(export.ndjson(recs))
        return 0
    problems = pack.verify(Path(a.root), a.release, packs=packs, scratch=scratch)
    _print(problems or [f"pack {a.release}: verified"])
    return 1 if problems else 0


DEFAULT_CPUS = 8  # `search run` without --cpus (and without a profile cap below it)


def _limits(a) -> runner.Limits:
    """The run's Limits (card W-10). Without a profile (--profile, else $BIOBOX_PROFILE): the flags
    alone, as before. With one: require_cap, nice, the thread cap for the current hour and the
    memory scope come from it, and --cpus / --mem-gb may lower them, never raise them (ValueError)."""
    cpus = DEFAULT_CPUS if a.cpus is None else a.cpus
    mem_gb = a.mem_gb
    extra: dict = {}
    path = a.profile or os.environ.get("BIOBOX_PROFILE")
    if path:
        from biobox import profile, runstep

        prof = profile.load(path)
        cap = runstep.choose_threads(prof, datetime.now(), runstep.read_cpu())  # card C-28: others' CPU, not the load
        if a.cpus is not None and a.cpus > cap:
            raise ValueError(f"--cpus {a.cpus} is above the profile's thread cap of {cap} for this hour")
        cpus = min(cpus, cap)
        mem_cap = runstep.size_bytes(prof.mem_scope)
        if mem_cap is not None:
            if mem_gb is not None and mem_gb * runner.GB > mem_cap:
                raise ValueError(f"--mem-gb {mem_gb:g} is above the profile's memory scope of {prof.mem_scope}")
            if mem_gb is None:
                mem_gb = mem_cap / runner.GB
        extra = {"require_cap": prof.require_mem_cap, "nice": prof.nice_monolithic}
    return runner.Limits(cpus=cpus, cpu_min=a.cpu_min, wall_s=a.wall_s, mem_gb=mem_gb, **extra)


def _pending_notice(run_id: str, info: Path) -> None:
    """`search run --wait` (card W-13): the run id, as soon as the run is pending, so a unit's log names it."""
    print(f"search run: {run_id} running; pending record {info}", file=sys.stderr, flush=True)


def _positive_prefs(controls_json: str | None) -> dict[bytes, tuple[str | None, str | None]]:
    """{protein_uid: (accession, header text)} for the positive controls of a run's `controls_json` (card W-20), so
    the FASTA's `example=` names the control's own contig: a label `ACCESSION:START-END` (or the GenBank form
    `ACCESSION:START..END`, card W-23) names its accession, any other label is header text (as `search hmm
    --positive` reads it). A control without a label names no contig."""
    out: dict[bytes, tuple[str | None, str | None]] = {}
    for p in (json.loads(controls_json) if controls_json else {}).get("positive", []):
        label = p.get("label")
        if label and p.get("uid"):
            try:
                span = parse_span(label)
            except SpanError:
                span = None
            out[bytes.fromhex(p["uid"])] = (span[0], None) if span else (None, label)
    return out


def _close_pending(final: Path) -> None:
    """Remove the pending record of a run that has closed (card W-20): `status` and `cancel` answer from the run's
    own record (`run.parquet`) from then on. Beside the run directory `final`, as `runner.pending_path` says."""
    runner.pending_path(Path(final).resolve()).unlink(missing_ok=True)


def _search(a) -> int:
    root = Path(a.root)
    if a.action == "run":
        unit = os.environ.get("BIOBOX_POD_UNIT")
        if unit and not a.wait:  # card W-13: decided by what `systemd-run --collect` allows (module docstring)
            print(f"refused: search run without --wait in the pod_run.sh unit {unit}: the unit ends when this "
                  "command returns, and systemd stops the detached runner with it. Add --wait: the run then keeps "
                  "the unit until it closes, and `bio-box search status` and `search cancel` work on it from "
                  "another shell", file=sys.stderr)
            return 2
        try:
            limits = _limits(a)
            tmp = _profile_tmp(a)
        except (ValueError, OSError) as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        # card W-03: the pack, the run, its pending record and its work directory, all outside the root
        packs = _packs(a)
        work = Path(a.work) if a.work else runner.default_work_base(tmp)
        out = Path(a.out) if a.out else None
        runs = None if out is not None else _runs(a)
        try:
            if not a.wait:
                print(runner.submit(root, Path(a.target), a.release, limits, work, a.fault, packs=packs, runs=runs,
                                    out=out))
                return 0
            row = runner.run(root, Path(a.target), a.release, limits=limits, work_base=work, argv=sys.argv,
                             run_id=bytes.fromhex(a.run_id) if a.run_id else None, fault=a.fault, packs=packs,
                             runs=runs, out=out, pending=not a.submitted, on_pending=_pending_notice)
        except RecipeError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        _close_pending(out if out is not None else runner.run_dir(runs, row["run_id"]))
        _print(export.ndjson([export.run_record(row)]) if a.json else [row["run_id"].hex() + "\t" + row["status"]])
        return 0 if row["status"] == "ok" else 1
    d = _run_path(a, a.target)
    if a.action == "status":
        try:
            print(json.dumps(runner.status(d)))
        except KeyError:
            print(f"refused: no run {a.target} (looked for {d} and its pending record)", file=sys.stderr)
            return 2
        return 0
    if a.action == "cancel":
        try:
            runner.cancel(d)
        except FileNotFoundError:
            if (d / "run.parquet").exists():  # card W-20: a closed run has no pending record; its own record answers
                print(f"run {a.target} is already closed ({runner.status(d)['status']}): nothing to cancel",
                      file=sys.stderr)
                return 0
            print(f"refused: no pending run {a.target} (looked for {runner.pending_path(d)})", file=sys.stderr)
            return 2
        return 0
    if not (d / "run.parquet").exists():
        print(f"refused: no closed run {a.target} (looked for {d})", file=sys.stderr)
        return 2
    try:
        tv = _run_tape_version(d, _packs(a))
    except RunPackMissing as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    hits = runner.read_hits(d, a.round)
    release = _run_release(d)
    manifest = export.release_manifest(root, release)
    if a.format:  # card W-05 (files: the manifest's, else the release_def stages: card W-12): one row per gene (tsv) or one record per protein (faa)
        if a.json:
            print("refused: --json and --format are exclusive", file=sys.stderr)
            return 2
        make = export.hits_tsv if a.format == "tsv" else export.hits_faa
        more = {"prefer": _positive_prefs(runner.read_run(d)["controls_json"])} if a.format == "faa" else {}
        try:
            _emit(make(root, hits, release=release, tape_version=tv, gene_call_version=a.gene_call_version,
                       manifest=manifest, library=a.library, threads=a.threads, **more), a.out,
                  "rows" if a.format == "tsv" else "sequences")
        except ValueError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        return 0
    occ = {}
    try:
        files = export.release_files(root, release, manifest)  # the run's release only (cards B-28, W-12)
    except (FileNotFoundError, ValueError) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    for o in nbr.occurrences(root, {h["target_uid"] for h in hits if h["target_kind"] == "protein"},
                             a.gene_call_version, tape_version=tv, files=files):
        occ.setdefault(o["protein_uid"], []).append(o)
    dims = export.contig_dims(root, release, {o["contig_uid"] for v in occ.values() for o in v}, manifest)
    recs = list(export.hit_records(hits, occ, a.top_n or None, a.members, tape_version=tv, dims=dims))
    if a.json:
        _print(export.ndjson(recs))
    else:
        cols = ["target_kind", "target_uid", "seq_bitscore", "bitscore", "evalue_zref", "env_start", "env_end",
                "n_occurrences", "curation_label"]
        _print(["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in recs])
        print(f"showing {len(recs)} of {len(hits)} rows" + (" (--top_n 0 shows all)" if len(recs) < len(hits) else ""),
              file=sys.stderr)
    return 0


def _anchor_run(a, source: str) -> Path | None:
    """The run directory `source` names (`_run_path`), when it holds a closed run; else None (protein uids)."""
    d = _run_path(a, source)
    return d if (d / "run.parquet").exists() else None


def _anchors(a, source: str, tape_version: str | None):
    """Anchor proteins from a recorded run's protein hits (a run directory, or a run id in the runs
    directory) or a comma-separated list of protein uids: (uids, tape version, the run's release or None),
    or a refusal message. A run's tape version is its pack's: a run whose pack the pack directory lacks is
    refused (card W-24)."""
    d = _anchor_run(a, source)
    if d is not None:
        try:
            run_tv = _run_tape_version(d, _packs(a))
        except RunPackMissing as e:
            return f"refused: {e}"
        uids = {h["target_uid"] for h in runner.read_hits(d) if h["target_kind"] == "protein"}
        if tape_version is not None and tape_version != run_tv:
            return f"refused: tape version {tape_version} is not the tape version {run_tv} that run {source} read"
        return uids, run_tv, _run_release(d)
    try:
        return {bytes.fromhex(x) for x in source.split(",") if x}, tape_version, None
    except ValueError:
        return f"refused: {source} is neither a closed run (in {_runs(a)} or as a directory) nor protein uids"


def _uid_file(path: str, flag: str) -> list[str] | str:
    """The protein uids of a file, one per line (blank lines and `#` lines skipped): `nbr --anchors-file` (card W-18)
    and `seq --proteins-file` (card W-23); a refusal message instead when it cannot be read, holds no uid, or holds a
    line that is not one."""
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as e:
        return f"refused: {flag} {path}: {e}"
    words = [x.strip() for x in lines if x.strip() and not x.lstrip().startswith("#")]
    bad = next((x for x in words if not re.fullmatch(r"[0-9a-fA-F]{32}", x)), None)
    if bad is not None or not words:
        return (f"refused: {flag} {path}: "
                + (f"{bad!r} is not a protein uid (32 hex digits)" if words else "no protein uids in it"))
    return words


def _nbr(a) -> int:
    root = Path(a.root)
    if a.action not in (None, "summarize", "flanks"):  # `nbr <run_id|uids>`: the first positional is the anchors
        if a.source is not None:
            print(f"refused: nbr: unknown action {a.action!r} (expected summarize or flanks)", file=sys.stderr)
            return 2
        a.action, a.source = None, a.action
    if a.anchors_file:  # card W-18: protein uids one per line, for a list too long for one argument
        if a.source or a.anchors:
            print("refused: --anchors-file is the anchors; give it without --anchors or a run id or uid list",
                  file=sys.stderr)
            return 2
        words = _uid_file(a.anchors_file, "--anchors-file")
        if isinstance(words, str):
            print(words, file=sys.stderr)
            return 2
        a.anchors = ",".join(words)
    if a.action in ("summarize", "flanks") and not (a.source or a.anchors):
        print(f"refused: nbr {a.action} needs a run id or a list of protein uids", file=sys.stderr)
        return 2
    if a.action is None and not (a.anchors or a.source):
        print("refused: nbr needs --anchors", file=sys.stderr)
        return 2
    source = a.source or a.anchors
    got = _anchors(a, source, a.tape_version)
    if isinstance(got, str):
        print(got, file=sys.stderr)
        return 2
    uids, tv, run_release = got
    if run_release is not None and a.release is not None and a.release != run_release:
        # card W-24: a run's anchors and tape version are its release's; another release's tapes, clusters and
        # contig_dim read with them gave an answer of neither release, with exit 0
        print(f"refused: run {source} ran against release {run_release}, not {a.release}: nbr reads a run with its "
              f"own release's files (drop --release, or name {run_release}). To place the run's proteins in release "
              f"{a.release}, give them as protein uids: write them one per line (bio-box sql, the last query of "
              f"AGENT_GUIDE.md section 5) and pass --anchors-file FILE --release {a.release}", file=sys.stderr)
        return 2
    release = a.release or run_release
    if release is None:
        print("refused: anchors given as protein uids need --release (the release whose files to read: its "
              "tapes, clusters and contig_dim)", file=sys.stderr)
        return 2
    try:
        files = export.release_files(root, release, export.release_manifest(root, release))  # card W-12
    except (FileNotFoundError, ValueError) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    if a.action == "flanks":  # card W-25
        return _nbr_flanks(a, root, uids, tv, release, files)
    policy = nbr.Policy.parse(a.window or "10g", a.strand, a.operon_gap, a.edge, a.circular)
    if a.action is None and a.format == "tsv":  # card W-05: one table per anchor locus
        if a.json:
            print("refused: --json and --format are exclusive", file=sys.stderr)
            return 2
        rows, _ = nbr.fetch(root, uids, policy, a.gene_call_version, tape_version=tv, files=files)
        bits = None
        if run_release is not None:
            bits = {}
            for h in runner.read_hits(_anchor_run(a, source)):
                if h["target_kind"] == "protein" and h.get("seq_bitscore") is not None:
                    bits[h["target_uid"]] = max(bits.get(h["target_uid"], h["seq_bitscore"]), h["seq_bitscore"])
        try:
            tables = export.locus_tables(root, rows, release=release, gene_call_version=a.gene_call_version,
                                         hit_bits=bits, manifest=files,
                                         library=a.library, threads=a.threads)
            if a.out is None:
                _print(line for k, (_, lines) in enumerate(tables) for line in ([""] if k else []) + lines)
            else:
                n_tables = n_rows = 0
                for name, lines in tables:
                    export.write_lines(lines, Path(a.out) / name)
                    n_tables, n_rows = n_tables + 1, n_rows + len(lines) - 3  # a header, the facts, the columns
                _wrote(a.out, f"{n_tables} tables, {n_rows} rows")  # card W-23
        except ValueError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        return 0
    if a.action == "summarize":
        by = [x.strip() for x in a.by.split(",") if x.strip()]
        if "pfam" in by and not any(True for _ in files.entries("l1_domains")):
            # card W-20: as `search hmm` does, say so once and leave the level out (it would report every neighbour
            # as unannotated, share 1.0)
            by = [x for x in by if x != "pfam"]
            if not by:
                print(f"refused: release {release} has no domain library, so --by pfam has nothing to summarize; "
                      "use --by nr50,nr90", file=sys.stderr)
                return 2
            print(f"note: release {release} has no domain library: the pfam level is left out; summarizing at "
                  f"{','.join(by)}", file=sys.stderr)
        try:
            s = nbr.summarize(root, uids, policy, by, release=release,
                              gene_call_version=a.gene_call_version, tape_version=tv, top_n=a.top,
                              manifest=files.manifest)
        except ValueError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        _emit(export.ndjson(export.nbr_summary_records(s)) if a.json else export.nbr_summary_tsv(s), a.out,
              "records" if a.json else "rows")
        return 0
    rows, counts = nbr.fetch(root, uids, policy, a.gene_call_version, tape_version=tv, files=files)
    if a.json:
        _print(export.ndjson(export.nbr_records(rows, a.top_n or None)))
    else:
        print(json.dumps(counts))
    return 0


def _nbr_flanks(a, root: Path, uids: set[bytes], tv: str | None, release: str, files) -> int:
    """`nbr flanks` (card W-25, `biobox.search.flanks`): the anchors' loci grouped by their flanks, written into
    `--out DIR` (new or empty, outside the root). 2: refused before anything is written."""
    from biobox.search import flanks

    def refuse(msg) -> int:
        print(f"refused: {msg}", file=sys.stderr)
        return 2

    try:
        window = flanks.parse_window(a.window or flanks.DEFAULT_WINDOW)
    except flanks.FlanksError as e:
        return refuse(e)
    level = a.level or flanks.DEFAULT_LEVEL
    group_by = a.group_by or flanks.DEFAULT_GROUP_BY
    distance = a.distance or "j"
    for flag, value, allowed in (("--level", level, flanks.LEVELS), ("--group-by", group_by, flanks.GROUP_BY),
                                 ("--distance", distance, flanks.DISTANCES), ("--export", a.export or "gff", ("gff",))):
        if value not in allowed:
            return refuse(f"{flag} {value}: expected one of {', '.join(allowed)}")
    cut = flanks.DEFAULT_CUT if a.cut is None else a.cut
    if not 0 <= cut <= 1:
        return refuse(f"--cut {cut}: a distance between 0 and 1")
    max_leaves = flanks.DEFAULT_MAX_LEAVES if a.max_loci is None else a.max_loci
    if max_leaves < 1:
        return refuse(f"--max-loci {max_leaves}: at least 1")
    if a.figure_rows is not None and not a.figure:
        return refuse("--figure-rows bounds the rows of the figure: give it with --figure")
    figure_rows = flanks.DEFAULT_FIGURE_ROWS if a.figure_rows is None else a.figure_rows
    if figure_rows < 1:
        return refuse(f"--figure-rows {figure_rows}: at least 1")
    if a.figure:  # card W-27
        why = flanks.figure_unavailable()
        if why:
            return refuse(why)
    if not a.out:
        return refuse("nbr flanks writes a directory of tables: give --out DIR (new or empty, outside the root)")
    out = Path(a.out)
    try:
        pack.outside_root(root, out, "the output directory")
    except ValueError as e:
        return refuse(e)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        return refuse(f"the output directory {out} is not empty: nbr flanks writes into a new or empty directory")
    try:
        counts, notes = flanks.run(root, uids, out, release=release, files=files, window=window, level=level,
                                   group_by=group_by, distance=distance, cut=cut, max_leaves=max_leaves,
                                   export_gff=a.export == "gff", library=a.library,
                                   gene_call_version=a.gene_call_version, tape_version=tv, threads=a.threads,
                                   tmp=_profile_tmp(a), figure=a.figure, figure_rows=figure_rows)
    except ValueError as e:
        return refuse(e)
    for line in notes:
        print(line, file=sys.stderr)
    _wrote(out / "groups.tsv", f"{counts['groups']} rows")
    _wrote(out / "context_types.tsv", f"{counts['context_rows']} rows")
    _wrote(out / "conservation.tsv", f"{counts['conservation_rows']} rows")
    _wrote(out / "groups", f"{counts['groups']} groups, {counts['leaves']} leaves, at most {counts['peak_leaves']} "
                           "in one tree")
    if a.figure:
        _wrote(out / "groups" / "*" / "figure.{svg,pdf}", f"{counts['figures']} figures, each as SVG and PDF")
    if a.export == "gff":
        _wrote(out / "gff", f"{counts['gff']} tables")
    _print([f"flanks\t{out}\t{counts['loci']} loci in {counts['groups']} groups (level {level}, window {window}g, "
            f"group by {group_by}, distance {distance}, cut {cut:g})"])
    return 0


POSITIVE_SHOWN = 20  # genes listed when a positive control matches no protein or several
_SPAN_0 = re.compile(r"^(\d+)-(\d+)$")  # START-END: 0-based, end exclusive (bio-box's own form)
_SPAN_GB = re.compile(r"^(\d+)\.\.(\d+)$")  # START..END: 1-based, both ends included (GenBank's form)
_SPAN_LIKE = re.compile(r"^[\d.\-]*\d[\d.\-]*$")  # digits, dots and dashes: meant as coordinates


class SpanError(ValueError):
    """Coordinates after `ACCESSION:` in neither form `parse_span` reads; the message shows both spellings."""


def parse_span(spec: str) -> tuple[str, int, int] | None:
    """(accession, start, end), 0-based and half-open, of a gene or region named `ACCESSION:START-END` (0-based, end
    exclusive, as bio-box prints) or `ACCESSION:START..END` (1-based, both ends included: GenBank's form, so
    `NC_001604:3171..5822` is `NC_001604:3170-5822`), card W-23. None when `spec` names no coordinates: no colon,
    whitespace in the accession, or text after the last colon other than digits with dots or dashes. SpanError for
    such text in neither form (`3171.-5822`, `3171...5822`, `3171`), and for a GenBank start of 0 or after its end."""
    acc, sep, tail = spec.rpartition(":")
    if not sep or not acc or re.search(r"\s", acc) or not _SPAN_LIKE.match(tail):
        return None
    m = _SPAN_0.match(tail)
    if m:
        return acc, int(m[1]), int(m[2])
    m = _SPAN_GB.match(tail)
    if m and 1 <= int(m[1]) <= int(m[2]):
        return acc, int(m[1]) - 1, int(m[2])
    nums = [int(x) for x in re.findall(r"\d+", tail)]
    a, b = (nums + nums)[:2] if len(nums) < 2 else nums[:2]
    a1 = max(a, 1)
    raise SpanError(f"{spec!r} is not a span bio-box reads; write one of the two forms. If {a} and {b} are GenBank's "
                    f"coordinates (1-based, both ends included): {acc}:{a1}..{max(b, a1)}, which is "
                    f"{acc}:{a1 - 1}-{max(b, a1)} in bio-box's own form. If they are bio-box's (0-based, end "
                    f"exclusive, as `find` and `search hits` print them): {acc}:{a}-{b}, which is GenBank "
                    f"{acc}:{a + 1}..{b}")


def span_forms(acc: str, start: int, end: int) -> tuple[str, str]:
    """A span in both forms `parse_span` reads: `ACC:START-END` (0-based, half-open) and `ACC:START..END` (GenBank)."""
    return f"{acc}:{start}-{end}", f"{acc}:{start + 1}..{end}"


def _release_files_or_refusal(a, root: Path):
    """(the release's files, its domain libraries) for `search hmm`: the files as `search hits` and `nbr` read them
    (`export.release_files`: the manifest, else the release_def stages; card W-12), the libraries those of its
    `l1_domains` files. ValueError, naming `bio-box release define`, for a release with neither."""
    files = export.release_files(root, a.release, export.release_manifest(root, a.release))
    return files, sorted({str(parts["library"]) for _, parts in files.entries("l1_domains")})


def _gene_line(g: dict) -> str:
    """One gene of a refusal's list: its span in bio-box's form with its strand, the GenBank form (card W-23), its
    protein uid and its contig's header."""
    zero, gb = span_forms(g["accession"] or g["contig_uid"].hex(), g["start"], g["end"])
    return f"  {zero}({g['strand']:+d})\t{gb}\t{g['protein_uid'].hex()}\t{g['header'] or ''}"


def _gene_lines(genes: list[dict]) -> list[str]:
    genes = sorted(genes, key=lambda g: (g["accession"] or "", g["contig_uid"], g["start"], g["end"], g["strand"]))
    lines = [_gene_line(g) for g in genes[:POSITIVE_SHOWN]]
    return lines + ([f"  ... and {len(genes) - POSITIVE_SHOWN} more genes"] if len(genes) > POSITIVE_SHOWN else [])


def _positive_uid(root: Path, files, spec: str, tape_version: str, gene_call_version: str,
                  threads: int | None = None) -> tuple[str | None, list[str]]:
    """The protein uid (hex) `--positive` names, and no lines; or None and the lines of a refusal that list what
    matched. `ACCESSION:START-END` (or the GenBank form `ACCESSION:START..END`, `parse_span`; card W-23) names the
    gene of a contig with that accession (`contig_label`, as `search hits --format tsv` shows it) and those
    coordinates; coordinates in neither form are refused with both spellings; any other text is a substring of
    contig headers (case-sensitive) and names every gene of those contigs. A gene is a CDS slot of the contig's tape
    in `tape_version` (the pack's).

    The contigs are found by a scan of `contig_label` (sorted by contig_uid: no table is sorted by accession or
    header, and a substring cannot prune); their tapes are read by contig_uid, which prunes (AGENTS §4.4)."""
    try:
        m = parse_span(spec)
    except SpanError as e:
        return None, [f"positive control: {e}"]
    genes, labels = find_genes(root, files, accession=m[0] if m else None, text=None if m else spec,
                               tape_version=tape_version, gene_call_version=gene_call_version, threads=threads)
    if m:
        exact = [g for g in genes if (g["start"], g["end"]) == (m[1], m[2])]
    else:
        exact = genes
    proteins = sorted({g["protein_uid"] for g in exact})
    if len(proteins) == 1:
        return proteins[0].hex(), []
    head = f"positive control {spec!r}"
    if proteins:
        return None, [f"{head} matched {len(proteins)} proteins at {len(exact)} genes of release {files.release}; "
                      "a positive control is one protein: give one of these genes as ACCESSION:START-END",
                      *_gene_lines(exact)]
    if not m:
        return None, [f"{head} matched no protein of release {files.release}: no contig header contains this text "
                      "(case-sensitive)"]
    if not labels:
        return None, [f"{head} matched no protein: no contig of release {files.release} has the accession "
                      f"{m[0]} (the accession column of `bio-box search hits --format tsv`)"]
    miss, near = _span_miss(genes, m)
    return None, [f"{head} matched no protein of release {files.release}: {miss}", *_gene_lines(near)]


def _span_miss(genes: list[dict], span: tuple[str, int, int]) -> tuple[str, list[dict]]:
    """Why no gene of the accession has the span (`parse_span`'s), and the genes that overlap it."""
    acc, start, end = span
    near = [g for g in genes if g["start"] < end and g["end"] > start]
    return (f"no gene of {acc} starts at {start} and ends at {end} (0-based, end exclusive, as `search hits "
            f"--format tsv` prints them; GenBank {start + 1}..{end})"
            + ("; the genes there:" if near else "; no gene overlaps that interval"), near)


def _contig_labels(root: Path, files, gene_call_version: str, threads: int | None, *, where: str | None = None,
                   contigs=()) -> dict[bytes, tuple[str, str]]:
    """{contig_uid: (accession, header)} as `search hits` shows them, for the contigs of `contig_label` that satisfy
    the SQL condition `where`, else for the contig uids `contigs`. A scan of `contig_label` when `where` is given
    (no table is sorted by accession or header); by uid it prunes (AGENTS §4.4). The scan is cheap at today's sizes
    (card W-15: 0.03 to 0.05 s, one row group per stage on `rN.3`, *measured*; about 1 s at 1 M label rows,
    *derived*; `planning/measurements/2026-11-12_positive_lookup_rN3.md`)."""
    import pyarrow as pa

    with export.open_store(root, release=files.release, gene_call_version=gene_call_version, manifest=files,
                           threads=threads, relations=("lab", "cidx", "dim")) as st:
        if where is not None:
            st.con.execute(f"CREATE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM lab WHERE {where}")
        else:
            st.con.register("_cin", pa.table({"contig_uid": pa.array(sorted(contigs), pa.binary(16))}))
            st.con.execute("CREATE TEMP TABLE _c AS SELECT DISTINCT contig_uid FROM _cin")
        export._annotate_contigs(st.con, with_ordinals=False)
        return {cu: (acc, header) for cu, acc, header in
                st.con.execute("SELECT contig_uid, accession, header FROM _lab").fetchall()}


def find_genes(root: Path, files, *, accession: str | None = None, text: str | None = None,
               protein: bytes | None = None, tape_version: str, gene_call_version: str,
               threads: int | None = None) -> tuple[list[dict], dict[bytes, tuple[str, str]]]:
    """The genes of one release that one of three keys names, with the label of every contig the key matched
    (`bio-box find` and `search hmm --positive`; card W-18). Give exactly one key: `accession` (every gene of the
    contigs of that accession: the CDS slots of their tapes in `tape_version`), `text` (the same for the contigs whose
    header contains it, case-sensitive), or `protein` (every gene that carries that protein uid, from
    `protein_occurrence`). Each gene is a dict with protein_uid, contig_uid, feature_uid, start, end (0-based,
    half-open), strand (+1 or -1), accession and header; sorted by accession, contig_uid, start, end, strand. The
    feature_uid is the gene's identity; no ordinal is returned (AGENTS §3.3). Returns (genes, labels), where labels
    is {contig_uid: (accession, header)} of the matched contigs, so a caller can tell a contig without a gene from no
    contig at all."""
    if sum(k is not None for k in (accession, text, protein)) != 1:
        raise ValueError("find_genes: give exactly one of accession, text, protein")
    genes: list[dict] = []
    if protein is not None:
        occ = {(o["contig_uid"], o["feature_uid"]): o
               for o in nbr.occurrences(root, {protein}, gene_call_version, files=files)}
        labels = _contig_labels(root, files, gene_call_version, threads, contigs={cu for cu, _ in occ})
        genes = [{"protein_uid": o["protein_uid"], "contig_uid": cu, "feature_uid": fu, "start": o["start"],
                  "end": o["end"], "strand": o["strand"], "accession": labels[cu][0], "header": labels[cu][1]}
                 for (cu, fu), o in occ.items() if cu in labels]
    else:
        where = (f"accession = {export._q(accession)}" if accession is not None
                 else f"contains(original_header, {export._q(text)})")
        labels = _contig_labels(root, files, gene_call_version, threads, where=where)
        for cu, tape in nbr._tapes(files, tape_version, set(labels)).items():
            for k, pu in enumerate(tape["protein_uids"]):
                if pu is not None:
                    genes.append({"protein_uid": pu, "contig_uid": cu, "feature_uid": tape["feature_uids"][k],
                                  "start": tape["start_bps"][k], "end": tape["end_bps"][k],
                                  "strand": tape["strands"][k], "accession": labels[cu][0], "header": labels[cu][1]})
    genes.sort(key=lambda g: (g["accession"] or "", g["contig_uid"], g["start"], g["end"], g["strand"]))
    return genes, labels


# card W-23 appended the gene's span in both forms `parse_span` reads: `span` (ACC:START-END, 0-based, half-open)
# and `genbank_span` (ACC:START..END, 1-based, both ends included)
FIND_COLUMNS = ("accession", "contig_uid", "start", "end", "strand", "feature_uid", "protein_uid", "span",
                "genbank_span")


def _find(a) -> int:
    """`bio-box find` (card W-18): the genes an accession (optionally with a span, in either form of `parse_span`) or a
    protein uid names, one TSV row each: accession, contig_uid, start, end (0-based, half-open), strand, feature_uid,
    protein_uid, and the span in both forms (card W-23). 2 when it names no gene or the span is in neither form."""
    root = Path(a.root)

    def refuse(msg, more=()) -> int:
        print(f"refused: {msg}", file=sys.stderr)
        for line in more:
            print(line, file=sys.stderr)
        return 2

    try:
        files, _ = _release_files_or_refusal(a, root)
    except (ValueError, OSError) as e:
        return refuse(e)
    tv = nbr.resolve_tape_version(a.gene_call_version, a.tape_version)
    if a.protein is not None:
        if not re.fullmatch(r"[0-9a-fA-F]{32}", a.protein):
            return refuse(f"--protein {a.protein!r} is not a protein uid (32 hex digits)")
        genes, _ = find_genes(root, files, protein=bytes.fromhex(a.protein), tape_version=tv,
                              gene_call_version=a.gene_call_version, threads=a.threads)
        if not genes:
            return refuse(f"no gene of release {files.release} carries protein {a.protein.lower()}")
    else:
        try:
            m = parse_span(a.accession)
        except SpanError as e:
            return refuse(f"--accession {e}")
        acc = m[0] if m else a.accession
        genes, labels = find_genes(root, files, accession=acc, tape_version=tv, gene_call_version=a.gene_call_version,
                                   threads=a.threads)
        if not labels:
            return refuse(f"no contig of release {files.release} has the accession {acc} (the accession column of "
                          "`bio-box search hits --format tsv`)")
        if m:
            miss, near = _span_miss(genes, m)
            genes = [g for g in genes if (g["start"], g["end"]) == (m[1], m[2])]
            if not genes:
                return refuse(miss, _gene_lines(near))
        elif not genes:
            return refuse(f"no gene of {acc} in tape version {tv} of release {files.release}")
    _print(["\t".join(FIND_COLUMNS)] + ["\t".join([g["accession"] or "", g["contig_uid"].hex(), str(g["start"]),
                                                   str(g["end"]), f"{g['strand']:+d}", g["feature_uid"].hex(),
                                                   g["protein_uid"].hex(),
                                                   *span_forms(g["accession"] or g["contig_uid"].hex(), g["start"],
                                                               g["end"])]) for g in genes])
    return 0


def _recipe_text(a, positive_uid: str) -> str:
    """The recipe `search hmm` runs and writes into DIR (recipe.load documents every key)."""
    seed = {"hmm_file": str(Path(a.hmm).resolve())}  # hashed by the file's sha256, not its path (card W-02)
    if a.model:
        seed["model"] = a.model
    name = re.sub(r"[^a-z0-9]+", "_", Path(a.hmm).stem.lower()).strip("_") or "seed"
    return yaml.safe_dump({
        "recipe": f"hmm_{name}", "written_against": a.release, "seed": [seed], "target": "unique",
        "tool": "hmmsearch", "iterations": {"max": 1, "inclusion_bits": float(a.bits)},
        "context": {"window": f"{a.neighbors}g"},
        "controls": {"positive": [{"protein_uid": positive_uid, "label": a.positive}]},
    }, sort_keys=False)


def _count(path: Path, prefix: str = "") -> int:
    with open(path, encoding="utf-8") as fh:
        return sum(1 for line in fh if line.startswith(prefix))


def _search_hmm(a) -> int:
    """`bio-box search hmm` (card W-07; module docstring): check everything that can be checked, run, then write DIR
    by calling the separate commands. 2: refused before the run; 1: the run ended other than `ok` (DIR then holds
    the recipe it ran); an export's own code when it fails (DIR is then not written; the run stays)."""
    root = Path(a.root)

    def refuse(msg, more=()) -> int:
        print(f"refused: {msg}", file=sys.stderr)
        for line in more:
            print(line, file=sys.stderr)
        return 2

    if not a.positive:
        return refuse("search hmm needs a positive control, a protein the search must find: --positive "
                      "ACCESSION:START-END (one gene: its accession, start and end as `bio-box search hits --format "
                      "tsv` prints them, 0-based, end exclusive), ACCESSION:START..END (the same gene in GenBank's "
                      "1-based form) or --positive 'TEXT' (a substring of the contig's header). A search without "
                      "one is refused (S-04)")
    if a.neighbors < 1:
        return refuse(f"--neighbors {a.neighbors}: give at least 1 gene on each side")
    out = Path(a.out)
    try:
        pack.outside_root(root, out, "the output directory")
    except ValueError as e:
        return refuse(e)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        return refuse(f"the output directory {out} is not empty: search hmm writes into a new or empty directory")
    try:
        limits = _limits(a)
        packs, runs = _packs(a), _runs(a)
        work = Path(a.work) if a.work else runner.default_work_base(_profile_tmp(a))
        tv = pack.manifest_tape_version(pack.read_manifest(pack.find(packs, a.release, root=root))[0])
        files, libraries = _release_files_or_refusal(a, root)
        export._Store(None, libraries).library(a.library)  # several libraries need --library, as `search hits` says
    except (ValueError, OSError) as e:  # a missing pack (pack.NoPack) names the `pack build` command
        return refuse(e)
    uid, lines = _positive_uid(root, files, a.positive, tv, a.gene_call_version, a.threads)
    if uid is None:
        return refuse(lines[0], lines[1:])
    text = _recipe_text(a, uid)
    try:  # in the foreground, with a pending record: `search status` and `search cancel` work on it (card W-13)
        row = runner.run(root, text, a.release, limits=limits, work_base=work, argv=sys.argv, packs=packs,
                         runs=runs, pending=True, on_pending=_pending_notice, gene_call_version=a.gene_call_version)
    except RecipeError as e:
        return refuse(e)
    _close_pending(runner.run_dir(runs, row["run_id"]))  # card W-20
    rid, status = row["run_id"].hex(), row["status"]
    levels = "pfam,nr50" if libraries else "nr50,nr90"
    window = f"{a.neighbors}g"
    tmp = out.parent / f".{out.name}.tmp-{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    global _QUIET
    quiet = _QUIET
    _QUIET = True  # card W-23: the commands below write into tmp; the lines name the files once DIR is in place
    try:
        (tmp / "recipe.yaml").write_text(text, encoding="utf-8")
        if status == "ok":
            (tmp / "README.txt").write_text(export.README_TEXT, encoding="utf-8")  # card W-20
            where = ["--root", str(root), "--gene-call-version", a.gene_call_version, "--runs", str(runs),
                     "--pack-dir", str(packs)]
            opts = (["--library", a.library] if a.library else []) + (["--threads", str(a.threads)] if a.threads else [])
            (tmp / "loci").mkdir()
            for argv in (["search", "hits", rid, *where, *opts, "--format", "tsv", "--out", str(tmp / "hits.tsv")],
                         ["search", "hits", rid, *where, *opts, "--format", "faa", "--out", str(tmp / "hits.faa")],
                         ["nbr", "summarize", rid, *where, "--window", window, "--by", levels, "--format", "tsv",
                          "--out", str(tmp / "neighbors.tsv")],
                         ["nbr", rid, *where, *opts, "--window", window, "--format", "tsv", "--out", str(tmp / "loci")]):
                code = main(argv)
                if code != 0:
                    print(f"search hmm: run {rid} closed ok, but `bio-box {' '.join(argv[:2])} ...` exited {code}; "
                          f"nothing was written to {out}; the run is {runner.run_dir(runs, rid)}", file=sys.stderr)
                    return code
            counts = (_count(tmp / "hits.tsv") - 1, _count(tmp / "hits.faa", ">"), len(list((tmp / "loci").iterdir())))
        os.replace(tmp, out)  # an empty DIR is replaced; a DIR that gained files meanwhile is not (OSError)
    finally:
        _QUIET = quiet
        shutil.rmtree(tmp, ignore_errors=True)
    if status == "ok":
        for name in ("hits.tsv", "hits.faa", "neighbors.tsv", "loci"):
            _wrote(out / name, _WRITTEN[str(tmp / name)])
    print(f"{rid}\t{status}")
    if status != "ok":
        print(f"search hmm: run {rid} ended {status}: {row.get('error')}; {out} holds the recipe it ran",
              file=sys.stderr)
        return 1
    _print([f"recipe\t{out / 'recipe.yaml'}", f"readme\t{out / 'README.txt'}", f"hits\t{out / 'hits.tsv'}\t{counts[0]} genes of {counts[1]} proteins",
            f"fasta\t{out / 'hits.faa'}", f"neighbors\t{out / 'neighbors.tsv'}\tlevels {levels}, window {window}",
            f"loci\t{out / 'loci'}\t{counts[2]} tables"])
    if not libraries:
        print(f"note\trelease {a.release} has no domain library: the neighbour summary is at the cluster levels "
              f"{levels} only, and pfam_architecture is -")
    return 0


def _hmm_parser() -> argparse.ArgumentParser:
    """`bio-box search hmm` (card W-07): its own parser, since `search`'s others take a recipe or a run."""
    h = argparse.ArgumentParser(prog="bio-box search hmm", description="search a release with an HMM file and write "
                                "the recipe, hits, FASTA and neighbour tables into --out (card W-07)")
    h.add_argument("--root", required=True)
    h.add_argument("--gene-call-version", default=GENE_CALL_VERSION)
    h.add_argument("--release", required=True)
    h.add_argument("--hmm", required=True, help="the HMM file (HMMER3); a file of several models needs --model")
    h.add_argument("--model", help="the model of --hmm to search with, by name or accession")
    h.add_argument("--bits", type=float, required=True, help="inclusion threshold in bits (iterations.inclusion_bits)")
    h.add_argument("--neighbors", type=int, default=5, help="genes on each side of a hit in the neighbour tables "
                                                            "(window <N>g; default 5)")
    h.add_argument("--positive", help="required: a protein the search must find, as ACCESSION:START-END (one gene, "
                                      "as `search hits --format tsv` prints it: 0-based, end exclusive), as "
                                      "ACCESSION:START..END (GenBank: 1-based, inclusive) or a substring of its "
                                      "contig's header")
    h.add_argument("--out", required=True, help="a new or empty directory outside the data root")
    h.add_argument("--profile", help="host profile YAML (default: $BIOBOX_PROFILE); caps cpus and memory, and gives "
                                     "the pack, runs and work directories")
    h.add_argument("--pack-dir", help=PACK_DIR_HELP)
    h.add_argument("--runs", help=RUNS_HELP)
    h.add_argument("--work", help="parent of the run's work directory (default: <profile tmp>/search/<user>/work "
                                  "with a profile, else $BIOBOX_SEARCH_WORK or /tmp/atd/search)")
    h.add_argument("--cpus", type=int, help=f"worker cpus (default {DEFAULT_CPUS}, or the profile's cap if lower)")
    h.add_argument("--mem-gb", type=float, help="memory cap in GB (default: the profile's mem_scope, else none)")
    h.add_argument("--cpu-min", type=float, default=300.0)
    h.add_argument("--wall-s", type=float, default=6 * 3600.0)
    h.add_argument("--library", help="domain library of pfam_architecture (default: the only one, else the pfam one)")
    h.add_argument("--threads", type=int, help=f"DuckDB threads of the exports (default {export.EXPORT_THREADS})")
    return h


def _export(a) -> int:
    """`export proteins` (card W-05): the release's distinct proteins as FASTA named by protein uid."""
    root = Path(a.root)
    manifest = export.release_manifest(root, a.release)
    if manifest is None:
        print(f"error: no manifests/release={a.release}.json under {root}; readers use the manifest's file list, "
              f"never a glob. Write it with: bio-box release manifest --root {root} --release {a.release}",
              file=sys.stderr)
        return 2
    _emit(export.protein_fasta(root, release=a.release, gene_call_version=a.gene_call_version, manifest=manifest,
                               threads=a.threads), a.out, "sequences")
    return 0


# --- card W-23: releases, seq, sql ---------------------------------------------------------------------------


def _level_rows(paths, level: str) -> int:
    """Rows at `level` of cluster membership files sorted by level first (`member_cluster`, `cluster_member`): from
    each row group's statistics on `level`; the column of a row group is read only where it holds several levels or
    has no statistics, so a file costs its footer and at most a row group per level boundary."""
    import pyarrow.compute as pc

    n = 0
    for p in paths:
        pf = pq.ParquetFile(p)
        k = pf.schema_arrow.get_field_index("level")
        for i in range(pf.metadata.num_row_groups):
            rg = pf.metadata.row_group(i)
            st = rg.column(k).statistics
            lo, hi = (st.min, st.max) if st is not None and st.has_min_max else (None, None)
            lo, hi = (x.decode() if isinstance(x, bytes) else x for x in (lo, hi))
            if lo is not None and lo == hi:
                n += rg.num_rows if lo == level else 0
            elif lo is not None and not lo <= level <= hi:
                continue
            else:
                n += pc.sum(pc.equal(pf.read_row_group(i, columns=["level"]).column(0), level)).as_py() or 0
    return n


def _distinct(paths, column: str) -> int:
    """Distinct values of one column of the files, read a row group at a time."""
    seen: set = set()
    for p in paths:
        pf = pq.ParquetFile(p)
        for i in range(pf.metadata.num_row_groups):
            seen.update(pf.read_row_group(i, columns=[column]).column(0).to_pylist())
    return len(seen)


def _release_genomes(files, release: str, dims: list) -> int | None:
    """The genomes of a release for `release_summary` (card C-21c); None without `contig_dim`."""
    if not dims:
        return None
    gdim = files.table("genome_dim", release=release)
    if gdim:
        return sum(pq.read_metadata(p).num_rows for p in gdim)
    member = files.table("genome_input")
    return _distinct(member or dims, "genome_uid")


def release_names(root: Path) -> list[str]:
    """Every release the root names: a manifest `manifests/release=<R>.json` or a release definition
    `release_def/release=<R>/`, so a release defined but not yet given a manifest is found too."""
    from biobox.registry.schemas import table

    names = set()
    man = Path(root) / "manifests"
    if man.is_dir():
        names |= {p.name[len("release="):-len(".json")] for p in man.glob("release=*.json")}
    defs = Path(root) / table("release_def").location
    if defs.is_dir():
        names |= {p.name[len("release="):] for p in defs.iterdir() if p.is_dir() and p.name.startswith("release=")}
    return sorted(names)


def release_summary(root: Path, release: str, packs: Path | None = None) -> dict:
    """What one release holds, for `bio-box releases` (card W-23), from its file list (the manifest, else its
    `release_def` stages: `catalog.ReleaseFiles`) and Parquet footers, never a scan of a large table: `contigs` the
    rows of its `contig_dim`, `proteins` its `nr90` members (`member_cluster`, else `cluster_member`, by row-group
    statistics on `level`), `genomes` every genome of the release (card C-21c): the rows of its `genome_dim` (one per
    genome, from footers), else, for a release built before card C-21b, the distinct `genome_uid` of its stages'
    `genome_input` (the catalog's `genome_member`; one column read: no footer holds a distinct count), else of its
    `contig_dim`, which names one genome per contig and so misses a genome whose every contig another genome also
    holds; None where the release has no `contig_dim`. `libraries` its domain libraries, `ncrna` its ncRNA
    producers, `lookups` whether it has the tables and column current readers use (cards B-27, B-39), and with a pack
    directory `pack`, the id of the release's pack there or None."""
    root = Path(root)
    rec: dict = {"release": release, "manifest": catalog.manifest_path(root, release).is_file()}
    try:
        files = catalog.ReleaseFiles.of(root, release)
    except (FileNotFoundError, ValueError, OSError) as e:
        return {**rec, "stages": None, "error": str(e)}
    dims = files.table("contig_dim", release=release)
    clus = files.table("member_cluster", release=release) or files.table("cluster_member", release=release)
    tok = files.table("tape_tokens", release=release)
    rec.update({
        "stages": list(files.stages),
        "genomes": _release_genomes(files, release, dims),
        "contigs": sum(pq.read_metadata(p).num_rows for p in dims) if dims else None,
        "proteins": _level_rows(clus, "nr90") if clus else None,
        "libraries": sorted({str(parts["library"]) for _, parts in files.entries("l1_domains")}),
        "ncrna": sorted({str(parts["producer"]) for _, parts in files.entries("l1_ncrna")}),
        "lookups": {"member_cluster": bool(files.table("member_cluster", release=release)),
                    "contig_label": bool(files.table("contig_label")),
                    "domain_family_tokens": bool(tok) and all(catalog.TAPE_DOMAIN_FAMILIES in pq.read_schema(p).names
                                                              for p in tok)},
    })
    if files.manifest is not None:  # card C-19: a library scanned on some of the release's stages only
        partial = {lib: missing for lib, _, missing in catalog.partial_libraries(files.manifest)}
        if partial:
            rec["partial_libraries"] = partial
    if packs is not None:
        try:
            rec["pack"] = pack.find(packs, release).name
        except pack.NoPack:
            rec["pack"] = None
        except ValueError:
            rec["pack"] = "several"
    return rec


def contig_slice(root: Path, files, contig_uid: bytes, start: int, end: int | None, strand: int,
                 circular: bool) -> tuple[bytes, int, dict]:
    """(the bases of `contig_uid`'s default variant from `start` to `end`, 0-based and half-open, on `strand`; the
    contig's length; what was read) for `bio-box seq --contig` (card W-23). Read by key through `SeqStore` over the
    release's files: one row group of `contig_index` and of `seq_locator` per key and the BGZF blocks of the slice
    (`reads`: the store's counters). `end` None is the contig's end; on a circular contig `end` may pass its length
    (across the origin). KeyError for a contig the release does not hold, ValueError for a slice outside it."""
    from biobox.l0seq.store import SeqStore

    store = SeqStore(root, files)
    try:
        n = store.length_of(store.default_digest(contig_uid))
        seq = store.get(contig_uid, start, n if end is None else end, strand, circular)
        reads = {"index_row_groups_read": store.index.row_groups_read,
                 "locator_row_groups_read": store.locator.row_groups_read, **store.counters()}
    finally:
        store.close()
    return seq, n, reads


def _seq(a) -> int:
    """`bio-box seq` (card W-23): proteins by uid, or one contig slice, as FASTA, read by key from the release's files
    (never by exporting the release). 2 when a key names nothing in the release."""
    root = Path(a.root)

    def refuse(msg) -> int:
        print(f"refused: {msg}", file=sys.stderr)
        return 2

    try:
        files = export.release_files(root, a.release, export.release_manifest(root, a.release))
    except (ValueError, OSError) as e:
        return refuse(e)
    if a.contig is None:
        if a.strand != 1:
            return refuse("--strand is for --contig; a protein has one strand")
        words = list(a.protein or [])
        if a.proteins_file:
            got = _uid_file(a.proteins_file, "--proteins-file")
            if isinstance(got, str):
                print(got, file=sys.stderr)
                return 2
            words = got
        bad = next((w for w in words if not re.fullmatch(r"[0-9a-fA-F]{32}", w)), None)
        if bad is not None:
            return refuse(f"--protein {bad!r} is not a protein uid (32 hex digits)")
        uids = list(dict.fromkeys(bytes.fromhex(w) for w in words))  # in the order given, each once
        seqs, _ = export.protein_seqs(root, uids, release=a.release, gene_call_version=a.gene_call_version,
                                      manifest=files, threads=a.threads)
        missing = [u.hex() for u in uids if u not in seqs]
        if missing:
            return refuse(f"{len(missing)} of {len(uids)} proteins are not in release {a.release}: "
                          f"{', '.join(missing[:20])}" + (" ..." if len(missing) > 20 else "") + "; nothing written")
        _emit((line for u in uids for line in (f">{u.hex()} release={a.release} length={seqs[u][0]}", seqs[u][1])),
              a.out, "sequences")
        return 0
    try:
        span = parse_span(a.contig)
    except SpanError as e:
        return refuse(f"--contig {e}")
    acc, start, end = span if span else (a.contig, 0, None)
    if re.fullmatch(r"[0-9a-fA-F]{32}", acc):  # a contig uid: accessions are never 32 hex digits
        cu = bytes.fromhex(acc)
        name = _contig_labels(root, files, a.gene_call_version, a.threads, contigs={cu}).get(cu, (None,))[0]
    else:  # by accession: a scan of contig_label, which is sorted by contig; cheap at today's sizes (card W-15)
        labels = _contig_labels(root, files, a.gene_call_version, a.threads, where=f"accession = {export._q(acc)}")
        if len(labels) != 1:
            return refuse(f"no contig of release {a.release} has the accession {acc}" if not labels else
                          f"{len(labels)} contigs of release {a.release} have the accession {acc}: "
                          f"{', '.join(sorted(c.hex() for c in labels))}; give one by its contig uid")
        cu, name = next(iter(labels)), acc
    circular = bool(export.contig_dims(root, a.release, {cu}, files).get(cu, {}).get("is_circular"))
    try:
        seq, n, _ = contig_slice(root, files, cu, start, end, a.strand, circular)
    except KeyError:
        return refuse(f"contig {cu.hex()} is not in release {a.release}")
    except ValueError as e:
        return refuse(f"--contig {a.contig}: {e} (0-based, end exclusive; an end past the contig's length only on a "
                      f"circular contig{'' if circular else ', and this one is linear'})")
    head = (f">{cu.hex()} release={a.release} accession={name or '-'} start={start} end={n if end is None else end} "
            f"strand={a.strand:+d} contig_len={n}")
    _emit([head, seq.decode("ascii")], a.out, "sequences")
    return 0


def sql_connection(root: Path, release: str, *, runs=(), tmp: Path | None = None, threads: int = 4,
                   memory_limit: str = "4GB"):
    """(a read-only connection to the release's catalog, its scratch directory) for `bio-box sql` (card W-23). The
    caller closes the connection and removes the scratch directory.

    The catalog is opened with `catalog.open_catalog`, so the sequence functions, `partners_of` and `token_of` work,
    from a copy in a scratch directory made under `tmp` (the profile's; default: the system's temporary directory;
    refused under the data root), which also holds DuckDB's temporary directory (card W-18: never under the root).
    The search runs in `runs` (run directories) are attached (`catalog.attach_search`). The protections it relies on,
    each enforced by DuckDB on this connection, not by reading the query's text:

    - the catalog is attached READ_ONLY and is the default database, so CREATE, INSERT, UPDATE, DELETE, DROP and
      ALTER on it fail;
    - external access is disabled (`enable_external_access = false`, which DuckDB refuses to turn back on while the
      database runs) with `allowed_paths` set to exactly the files the catalog's views and the runs read (the
      manifest's file list, `catalog_files`, and each run's Parquet files), so no other file on the host can be read
      or written: COPY ... TO, EXPORT, ATTACH, INSTALL and LOAD of any path fail, and `allowed_paths` and
      `allowed_directories` cannot be changed any more. DuckDB writes even a listed path through a temporary file
      beside it, which is not allowed, so a listed file cannot be overwritten either (both measured on wilkcomp02,
      DuckDB 2.0.0, 2026-11-12, and tested). The one directory DuckDB still writes is its own temporary directory, in
      the scratch directory, outside the root;
    - extensions are neither installed nor loaded automatically.

    `bio-box sql` also runs only statements that DuckDB's own parser classifies as SELECT (`_sql`)."""
    root = Path(root)
    path = catalog.catalog_path(root, release)
    if not path.is_file():
        raise FileNotFoundError(f"release {release} has no catalog at {path}; it is built by: bio-box release catalog "
                                f"--root {root} --release {release}")
    if tmp is not None:
        Path(tmp).mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="biobox-sql-", dir=tmp))
    try:
        pack.outside_root(root, scratch, "the scratch directory of bio-box sql")
        (scratch / "duckdb" / "profile").mkdir(parents=True)
        con = catalog.open_catalog(path, local_dir=scratch / "catalog", threads=threads, memory_limit=memory_limit,
                                   temp_directory=scratch / "duckdb", runs=list(runs))
    except BaseException:
        shutil.rmtree(scratch, ignore_errors=True)
        raise
    try:
        built = Path(con.execute("SELECT value FROM catalog_info WHERE key = 'root'").fetchall()[0][0])
        allowed = {str(built / p) for (p,) in con.execute("SELECT path FROM catalog_files").fetchall()}
        for d in runs:
            for name in catalog.SEARCH_RUN_FILES.values():
                if (Path(d) / name).is_file():
                    allowed |= {str(Path(d) / name), str((Path(d) / name).resolve())}
        con.execute("SET autoinstall_known_extensions = false")
        con.execute("SET autoload_known_extensions = false")
        if allowed:
            con.execute(f"SET allowed_paths = [{', '.join(catalog._sql_str(p) for p in sorted(allowed))}]")
        con.execute("SET enable_external_access = false")
    except BaseException:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
        raise
    return con, scratch


def _sql_runs(a, tmp: Path | None) -> tuple[list[Path], int] | str:
    """The closed runs of the release `bio-box sql` attaches: each `--runs` value is a run directory or a runs
    directory (its `run=*` directories); with neither, the runs directory of the host profile, when there is one.
    Returns (runs of the release, runs of other releases left out), or a refusal message."""
    dirs = [Path(x) for x in a.runs] if a.runs else ([runner.default_runs_dir(tmp)] if tmp is not None else [])
    found: list[Path] = []
    for d in dirs:
        if (d / "run.parquet").is_file():
            found.append(d)
        elif d.is_dir():
            found += sorted(p.parent for p in d.glob("run=*/run.parquet"))
        elif a.runs:
            return f"refused: --runs {d}: neither a run directory nor a runs directory"
    mine = [d for d in found if _run_release(d) == a.release]
    return mine, len(found) - len(mine)


def _sql_cell(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (bytes, bytearray)):
        return bytes(v).hex()
    return str(v).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def _scan_name(path: str, built: Path) -> str:
    """The registry table a scanned file belongs to (`catalog.classify`), or the search run file's table."""
    p = Path(path)
    try:
        hit = catalog.classify(p.relative_to(built).as_posix())
        if hit is not None:
            return hit[0].name
    except (ValueError, catalog.LayoutError):
        pass
    return {f: t for t, f in catalog.SEARCH_RUN_FILES.items()}.get(p.name, p.name)


def _sql(a) -> int:
    """`bio-box sql` (card W-23): SELECT statements on the release's catalog (`sql_connection`), each statement's rows
    as a TSV table (a blank line between tables) or, for one statement, JSON (one object per row), to stdout or
    `--out` (refused under the data root); on stderr, per statement, its rows and the row groups it read per table
    (`catalog.profile_query`). 2: refused before anything ran; 1: a statement failed."""
    root = Path(a.root)

    def refuse(msg) -> int:
        print(f"refused: {msg}", file=sys.stderr)
        return 2

    if (a.file is None) == (a.query is None):
        return refuse("bio-box sql takes a query as a FILE or with -e QUERY, one of the two")
    try:
        text = Path(a.file).read_text(encoding="utf-8") if a.file is not None else a.query
    except (OSError, UnicodeDecodeError) as e:
        return refuse(f"{a.file}: {e}")
    try:
        stmts = duckdb.extract_statements(text)
    except duckdb.Error as e:
        return refuse(f"the query does not parse: {e}")
    if not stmts:
        return refuse("no statement in the query")
    for k, s in enumerate(stmts, 1):  # by DuckDB's parser, before anything runs
        kind = str(s.type).rsplit(".", 1)[-1]
        if kind != "SELECT":
            return refuse(f"statement {k} is a {kind} statement ({' '.join(s.query.split())[:80]}); bio-box sql runs "
                          "SELECT statements only (also SHOW, DESCRIBE, SUMMARIZE and FROM ...), on a read-only "
                          "connection: nothing is written, attached, installed, loaded or set")
    if a.format == "json" and len(stmts) > 1:
        return refuse(f"--format json takes one statement; the query has {len(stmts)}")
    try:
        if a.out is not None:
            pack.outside_root(root, Path(a.out), "--out")
        tmp = _profile_tmp(a)
        got = _sql_runs(a, tmp)
        if isinstance(got, str):
            print(got, file=sys.stderr)
            return 2
        runs, others = got
        con, scratch = sql_connection(root, a.release, runs=runs, tmp=tmp, threads=a.threads,
                                      memory_limit=a.memory_limit)
    except (ValueError, OSError) as e:
        return refuse(e)
    if others:
        print(f"note: {others} runs of other releases left out", file=sys.stderr)
    try:
        built = Path(con.execute("SELECT value FROM catalog_info WHERE key = 'root'").fetchall()[0][0])
        out: list[str] = []
        for k, s in enumerate(stmts, 1):
            try:
                p = catalog.profile_query(con, s.query, profile_dir=scratch / "duckdb" / "profile")
            except duckdb.Error as e:
                print(f"sql: statement {k} failed: {e}", file=sys.stderr)
                return 1
            if a.format == "json":
                out += [json.dumps({c: export._plain(v) for c, v in zip(p["columns"], row)}, default=str)
                        for row in p["result"]]
            else:
                out += ([""] if k > 1 else []) + ["\t".join(p["columns"])]
                out += ["\t".join(_sql_cell(v) for v in row) for row in p["result"]]
            reads = [f"{_scan_name(sc['files'][0], built)} {sc['row_groups_scanned']} of "
                     f"{'?' if sc['row_groups_total'] is None else sc['row_groups_total']} row groups"
                     + (f" ({sc['lookup']})" if sc.get("lookup") else "") for sc in p["scans"] if sc["files"]]
            seq = p["bgzf"]["blocks_read"], sum(p["index"].values())
            print(f"sql: statement {k}: {len(p['result'])} rows; read "
                  + (", ".join(reads) or "no Parquet row group")
                  + (f"; sequence: {seq[0]} BGZF blocks, {seq[1]} index row groups" if any(seq) else ""),
                  file=sys.stderr)
    finally:
        con.close()
        shutil.rmtree(scratch, ignore_errors=True)
    _emit(out, a.out, "records" if a.format == "json" else "rows")
    return 0


RELEASES_COLUMNS = ("release", "manifest", "stages", "genomes", "contigs", "proteins", "libraries", "ncrna",
                    "member_cluster", "contig_label", "domain_family_tokens")


def _releases(a) -> int:
    """`bio-box releases` (card W-23): one line per release of the root (`release_summary`); `--json`, one object per
    line. With a profile or `--pack-dir`, the pack each release has there."""
    root = Path(a.root)
    packs = _packs(a) if (a.pack_dir or a.profile or os.environ.get("BIOBOX_PROFILE")) else None
    recs = [release_summary(root, r, packs) for r in release_names(root)]
    if not recs:
        print(f"note: no release under {root} (no manifests/release=*.json and no release_def/release=*)",
              file=sys.stderr)
    for r in recs:
        if "error" in r:
            print(f"note: release {r['release']}: {r['error']}", file=sys.stderr)
    if a.json:
        _print(json.dumps(r, sort_keys=False) for r in recs)
        return 0

    def text(v) -> str:
        if v is None:
            return "-"
        if isinstance(v, bool):
            return "yes" if v else "no"
        if isinstance(v, list):
            return ",".join(v) or "-"
        return str(v)

    cols = RELEASES_COLUMNS + (("pack",) if packs is not None else ())
    rows = []
    for r in recs:
        flat = {**r, **r.get("lookups", {})}
        libs = [f"{x}(not {','.join(r['partial_libraries'][x])})" if x in r.get("partial_libraries", {}) else x
                for x in r.get("libraries", [])]
        rows.append("\t".join(text(libs if c == "libraries" else flat.get(c)) for c in cols))
    _print(["\t".join(cols), *rows])
    return 0


def _run_step(a) -> int:
    from biobox import profile, runstep

    return runstep.run_step(profile.load(a.profile), Path(a.root), a.stage, a.step, a.command,
                            nice=a.nice, threads=a.threads)


RELEASE_STEPS = ("define", "cluster", "tokens", "postings", "dim", "manifest", "catalog", "verify")
CLUSTER_ACTIONS = ("plan", "step", "finish", "redo")  # `release cluster plan|step|finish|redo` (card B-41)
STAGE_STEPS = ("inputs", "ingest", "genecall", "igr", "tapes", "domains", "ncrna")
DOMAIN_ACTIONS = ("plan", "shard", "merge")  # `stage domains plan|shard|merge` (card B-24)
INGEST_ACTIONS = ("pass1", "plan", "batch", "merge")  # `stage ingest pass1|plan|batch|merge` (card B-30)
NCRNA_ACTIONS = ("plan", "shard", "merge")  # `stage ncrna plan|shard|merge` (card E-07)
STEP_ACTIONS = {"domains": DOMAIN_ACTIONS, "ingest": INGEST_ACTIONS, "ncrna": NCRNA_ACTIONS}
PACK_DIR_HELP = ("the pack directory, outside the data root (card W-03); default: <profile tmp>/packs, the host's one "
                 "pack directory (card W-22), with a host profile, else packs beside $BIOBOX_SEARCH_WORK or "
                 "/tmp/atd/search")
RUNS_HELP = ("the runs directory, outside the data root (card W-03); default: <profile tmp>/search/<user>/runs with "
             "a host profile, else runs in $BIOBOX_SEARCH_WORK or /tmp/atd/search")
PROFILE_DIRS_HELP = "host profile YAML (default: $BIOBOX_PROFILE); its tmp holds the default pack and runs directories"
OUT_HELP = ("output root (card C-08): every output goes under it at its root-relative path; default "
            "$BIOBOX_STAGING when set (bio-box run-step sets it and publishes from it), else --root")


def _csv(s: str | None) -> list[str] | None:
    return None if s is None else [x.strip() for x in s.split(",")]


def _out(a, root: Path) -> Path:
    """Where a builder step writes (card C-08): `--out`, else `$BIOBOX_STAGING` (set by run-step, which
    publishes from it), else the root itself."""
    if a.out:
        return Path(a.out)
    staging = os.environ.get("BIOBOX_STAGING")
    return Path(staging) if staging else root


def _release_step(a, root: Path, stages: list[str], out: Path) -> dict:
    tmp = Path(a.tmp) if a.tmp else None
    if a.step == "cluster":
        from biobox.l2 import cluster

        if a.action == "plan":  # clustering in parts (card B-41)
            if tmp is None:
                raise ValueError("release cluster plan needs --tmp (the scratch directory of the work)")
            return cluster.plan(root, a.release, a.gene_call_version, stages, tmp, parts=a.parts,
                                part_cpu_min=a.part_cpu_min, one_process_cpu_min_max=a.one_process_cpu_min_max,
                                out=out)
        return cluster.run(root, a.release, a.gene_call_version, stages, threads=a.threads, tmp=tmp, out=out)
    if a.step == "tokens":
        from biobox.l3 import tokens

        return tokens.run(root, a.release, stages, a.gene_call_version, libraries=_csv(a.libraries),
                          producers=_csv(a.producers), tmp=tmp, out=out)
    if a.step == "postings":
        from biobox.l3 import postings

        return postings.run(root, a.release, stages, a.gene_call_version, producers=_csv(a.producers), tmp=tmp,
                            out=out)
    if a.step == "dim":
        from biobox import dim

        return dim.run(root, a.release, stages, tmp=tmp, out=out)
    from biobox import catalog

    if a.step == "catalog":
        path = catalog.manifest_path(root, a.release)
        if not path.exists():
            raise FileNotFoundError(f"{path}: no manifest of release {a.release}; the catalog is built from it. "
                                    f"Write it with: bio-box release manifest --root {root} --release {a.release}")
        man = catalog.read_manifest(path)
        got = catalog.generate_catalog(man, root=root, out=out, tape_version=a.tape_version)
        return {"catalog": got.name, "manifest_hash": catalog.manifest_hash(man)}
    path = catalog.manifest_path(root, a.release)
    same_root = out.resolve() == root.resolve()
    if not same_root:  # refused before the manifest is built (write_manifest refuses again)
        for p in (path, catalog.manifest_path(out, a.release)):
            if p.exists():
                raise FileExistsError(f"{p}: the manifest of release {a.release} is written once; "
                                      "refusing to overwrite")
    man = catalog.build_manifest(root, a.release, stages)
    res = {"files": len(man["files"]), "artifacts": len(man["artifacts"]), "manifest_hash": catalog.manifest_hash(man)}
    if same_root and path.exists():  # into the root itself, as before card C-08: the same manifest is a no-op
        if catalog.manifest_hash(catalog.read_manifest(path)) != catalog.manifest_hash(man):
            raise ValueError(f"{path}: a different manifest of release {a.release} exists; refusing to overwrite")
        return res
    catalog.write_manifest(man, out=out)
    return res


def _verify(a) -> int:
    """`release verify` (card C-13): the release's manifest against the root (`catalog.verify_manifest`: every listed
    file there with its size and sha256, and no unlisted file in the release's scope) and, with `--previous P`, P's
    manifest against the release's (`contracts.check_append_only`, R-09: every file of P outside P's own `release=`
    partition is listed with the same sha256). Writes nothing. Prints each problem, then one line per check; exit 0
    when there is no problem, else 1 (a missing manifest is a problem). Warnings, which leave the exit code as it is,
    go to stderr: a lookup table the release lacks (card W-17), and a domain library the release has files of for
    some of its stages only (card C-19, `catalog.partial_libraries`)."""
    from biobox import contracts

    root = Path(a.root)

    def load(release: str) -> dict | None:
        path = catalog.manifest_path(root, release)
        if path.is_file():
            return catalog.read_manifest(path)
        print(f"release {release}: no manifest at {path} (bio-box release manifest --root {root} --release {release} "
              f"writes it)")
        return None

    man = load(a.release)
    if man is None:
        return 1
    problems = catalog.verify_manifest(man, root)
    _print(problems)
    print(f"release {a.release}: {'verified' if not problems else f'{len(problems)} problems'} (stages "
          f"{','.join(man['stages'])}, {len(man['files'])} files, {len(man.get('artifacts', []))} artifacts, manifest "
          f"hash {catalog.manifest_hash(man)})")
    print(catalog.describe_built_by(catalog.built_by(root, a.release, man["stages"])))  # card C-23: information only
    failed = bool(problems)
    have = {f["table"] for f in man["files"]}
    for lookup, source in (("member_cluster", "cluster_member"), ("contig_label", "contig_occurrence")):
        if source in have and lookup not in have:  # a release built before card B-27 (card W-17): a warning only
            print(f"warning: release {a.release} has no {lookup} files (built before card B-27); readers fall back to "
                  f"{source}, which cannot prune. A release built with the current code has the table "
                  f"(card C-17 builds one over these stages)", file=sys.stderr)
    for lib, scanned, unscanned in catalog.partial_libraries(man):  # card C-19: a warning only
        print(f"warning: release {a.release} has domain files of library {lib} for stages {','.join(scanned)} but "
              f"not {','.join(unscanned)} (not scanned when the release was built); its {lib} family tokens cover the "
              f"genes of {','.join(scanned)} only. A release built after those stages are scanned covers them all",
              file=sys.stderr)
    if a.previous:
        prev = load(a.previous)
        if prev is None:
            return 1
        rewritten = contracts.check_append_only(prev, man)
        _print(rewritten)
        print(f"append-only against {a.previous}: {'ok' if not rewritten else f'{len(rewritten)} problems'}")
        failed = failed or bool(rewritten)
    return 1 if failed else 0


def _release(a) -> int:
    """One release-level step over the release's stages (card B-25). `define` writes `release_def`
    once (each stage must have been ingested); a builder step given `--stages` defines the release if
    it is not yet defined; without `--stages` the stages come from `release_def`. `manifest` and
    `catalog` need the release defined; `catalog` needs its manifest. Inputs are read from `--root`,
    outputs written under `_out` (card C-08). `verify` writes nothing (`_verify`)."""
    from biobox import catalog, io

    if a.step == "verify":
        return _verify(a)
    root = Path(a.root)
    if a.action is not None and a.step != "cluster":
        print(f"refused: release {a.step} has no action {a.action} (only cluster: {', '.join(CLUSTER_ACTIONS)})",
              file=sys.stderr)
        return 2
    if a.step == "cluster" and a.action in ("step", "finish", "redo"):
        return _cluster_action(a, root)
    try:
        out = _out(a, root)
        given = _csv(a.stages)
        if a.step == "define":
            if given is None:
                raise ValueError("release define needs --stages")
            stages = catalog.check_stages(given)
            missing = [s for s in stages if not io.part_files(root / f"l0seq/contig_index/stage={s}")]
            if missing:
                raise FileNotFoundError(f"stages {missing} have no contig_index under {root}: not ingested")
            res: dict = {}
            stages = catalog.define_release(root, a.release, stages, out=out)
        else:
            if given is None or a.step in ("manifest", "catalog"):
                stages = catalog.release_stages(root, a.release)
                if given is not None and catalog.check_stages(given) != stages:
                    raise ValueError(f"release_def of {a.release} lists stages {stages}, not {given}")
            else:
                stages = catalog.verify_release_def(root, a.release, given)  # the builder defines it
            res = _release_step(a, root, stages, out)
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    head = {} if a.action is None else {"action": a.action}
    print(json.dumps({"release": a.release, "stages": stages, "step": a.step, **head, "result": res}, sort_keys=True,
                     default=str))
    return 0


def _cluster_action(a, root: Path) -> int:
    """`release cluster step|finish|redo` (card B-41): one step of clustering in parts after its plan. The stages
    are the plan's; `--tmp` is the scratch directory the work directory lives in (the same for every step)."""
    from biobox.l2 import cluster

    try:
        if not a.tmp:
            raise ValueError(f"release cluster {a.action} needs --tmp (the scratch directory of the plan's work)")
        if a.action != "finish" and not a.name:
            raise ValueError(f"release cluster {a.action} needs --name (a step of {cluster.PLAN_STEP}'s sequence)")
        tmp = Path(a.tmp)
        if a.action == "redo":  # writes nothing into the root or staging: not a run-step
            res = cluster.redo(root, a.release, a.name, tmp, a.threads)
        else:
            out = _out(a, root)
            if a.action == "step":
                res = cluster.step(root, a.release, a.name, tmp, a.threads, out=out)
            else:
                res = cluster.finish(root, a.release, tmp, a.threads, out=out)
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    print(json.dumps({"release": a.release, "step": "cluster", "action": a.action, "result": res}, sort_keys=True,
                     default=str))
    return 0


def _json_list(path: str) -> list[dict]:
    data = json.loads(Path(path).read_text())
    if not isinstance(data, list) or not all(isinstance(x, dict) for x in data):
        raise ValueError(f"{path}: expected a JSON list of objects")
    return data


def _genomes(path: str):
    """`--genomes FILE`: a JSON list of `GenomeInput` fields, contig uids as hex."""
    from biobox.l0seq.inputs import GenomeInput

    out = []
    for g in _json_list(path):
        kw = dict(g)
        kw["contig_uids"] = tuple(bytes.fromhex(u) for u in g["contig_uids"])
        kw["circular"] = tuple(bytes.fromhex(u) for u in g.get("circular", ()))
        out.append(GenomeInput(**kw))
    return out


def _source_objects(path: str):
    """`--objects FILE`: a JSON list of `SourceObject` fields (`retrieved_at` as ISO 8601); a relative
    `path` is taken from the directory of FILE."""
    from biobox.l0seq import SourceObject

    out = []
    for o in _json_list(path):
        kw = dict(o)
        p = Path(kw["path"])
        kw["path"] = p if p.is_absolute() else Path(path).parent / p
        if "retrieved_at" in kw:
            kw["retrieved_at"] = datetime.fromisoformat(kw["retrieved_at"])
        out.append(SourceObject(**kw))
    return out


def _stage_step(a, root: Path, out: Path) -> dict:
    tmp = Path(a.tmp) if a.tmp else None
    mfb = {} if a.max_file_bytes is None else {"max_file_bytes": a.max_file_bytes}
    if a.action is not None and a.action not in STEP_ACTIONS.get(a.step, ()):
        raise ValueError(f"stage {a.step}: {a.action} is not one of its actions "
                         f"({'; '.join(f'{s}: ' + ', '.join(v) for s, v in STEP_ACTIONS.items())})")
    if a.step == "inputs":
        from biobox.l0seq.inputs import write_genome_input

        if not a.genomes:
            raise ValueError("stage inputs needs --genomes")
        return write_genome_input(root, a.stage, _genomes(a.genomes), out=out, **mfb)
    if a.step == "ingest":
        from biobox.l0seq import ingest

        if not a.objects:
            raise ValueError("stage ingest needs --objects")
        kw = {k: v for k, v in (("shard_bases", a.shard_bases), ("line_width", a.line_width)) if v is not None}
        if a.action is None:
            return ingest(_source_objects(a.objects), root, a.stage, tmp=tmp, out=out, **kw, **mfb)
        from biobox.l0seq import store  # in batch processes (card B-30)

        need = {"pass1": ("shard", "pass1_bytes"), "batch": ("batch",)}.get(a.action, ())
        missing = [f"--{n.replace('_', '-')}" for n in need if getattr(a, n) is None]
        if missing:
            raise ValueError(f"stage ingest {a.action} needs {', '.join(missing)}")
        objs = _source_objects(a.objects)
        shared = {"runs_tmp": Path(a.runs_tmp) if a.runs_tmp else None}
        if a.action == "pass1":
            return store.ingest_pass1(objs, root, a.stage, a.shard, a.pass1_bytes, tmp=tmp, out=out, **kw, **shared)
        if a.action == "plan":
            return store.ingest_plan(objs, root, a.stage, tmp=tmp, out=out, **kw, **shared)
        if a.action == "batch":
            return store.ingest_batch(objs, root, a.stage, a.batch, tmp=tmp, out=out, **shared)
        return store.ingest_merge(objs, root, a.stage, tmp=tmp, out=out, **mfb, **shared)
    if a.step == "genecall":
        from biobox.l1.genecall import run

        return run(root, a.stage, version=a.gene_call_version, prefix_width=a.prefix_width, tmp=tmp,
                   workers=a.workers, seq_stages=_csv(a.seq_stages), out=out, **mfb)
    if a.step == "igr":
        from biobox.l1 import igr

        return igr.run(root, a.stage, gene_call_version=a.gene_call_version, producers=_csv(a.producers), tmp=tmp,
                       out=out, **mfb)
    if a.step == "tapes":
        from biobox.l3 import tapes

        kw = {} if a.max_slots is None else {"max_slots": a.max_slots}
        return tapes.run(root, a.stage, gene_call_version=a.gene_call_version, producers=_csv(a.producers), tmp=tmp,
                         out=out, **kw, **mfb)
    if a.step == "ncrna":  # cards E-01, E-07
        from biobox.l1 import ncrna

        need = {None: ("cm",), "plan": ("cm", "scan_bases"), "shard": ("cm", "shard"), "merge": ("cm",)}[a.action]
        missing = [f"--{n.replace('_', '-')}" for n in need if getattr(a, n) in (None, "")]  # shard 0 is given
        if missing:
            raise ValueError(f"stage ncrna{'' if a.action is None else ' ' + a.action} needs {', '.join(missing)}")
        kw = {"model_set": a.ncrna_set}
        cm = Path(a.cm)
        if a.action == "plan":
            return ncrna.plan(root, a.stage, cm, a.scan_bases, gene_call_version=a.gene_call_version, tmp=tmp,
                              out=out, **kw)
        if a.action == "shard":
            return ncrna.scan_shard(root, a.stage, cm, a.shard, a.cpus, seq_stages=_csv(a.seq_stages), tmp=tmp,
                                    out=out, **kw)
        if a.action == "merge":
            return ncrna.merge(root, a.stage, cm, tmp=tmp, out=out, **kw, **mfb)
        return ncrna.run(root, a.stage, cm, a.cpus, gene_call_version=a.gene_call_version, tmp=tmp, out=out,
                         seq_stages=_csv(a.seq_stages), **kw, **mfb)
    from biobox.l1 import domains

    need = {None: ("library", "library_version", "hmm", "cutoff"),  # the whole scan in this process
            "plan": ("library", "library_version", "hmm", "cutoff", "shard_residues"),  # card B-24
            "shard": ("library", "library_version", "hmm", "shard"),
            "merge": ("library", "library_version")}[a.action]
    missing = [f"--{n.replace('_', '-')}" for n in need if getattr(a, n) in (None, "")]  # shard 0 is given
    if missing:
        raise ValueError(f"stage domains{'' if a.action is None else ' ' + a.action} needs {', '.join(missing)}")
    if a.action == "plan":
        return domains.plan(root, a.stage, a.gene_call_version, a.library, a.library_version, Path(a.hmm), a.cutoff,
                            a.shard_residues, out=out)
    if a.action == "shard":
        return domains.scan_shard(root, a.stage, a.library, a.library_version, Path(a.hmm), a.shard, a.cpus,
                                  tmp=tmp, out=out)
    if a.action == "merge":
        return domains.merge(root, a.stage, a.library, a.library_version, tmp=tmp, out=out, **mfb)
    return domains.run(root, a.stage, a.gene_call_version, a.library, a.library_version, Path(a.hmm), a.cutoff,
                       cpus=a.cpus, out=out, **mfb)


def _stage(a) -> int:
    """One per-stage builder step (card C-08): inputs are read from `--root`, outputs written under
    `--out`, else `$BIOBOX_STAGING`, else the root."""
    root = Path(a.root)
    try:
        res = _stage_step(a, root, _out(a, root))
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    head = {"stage": a.stage, "step": a.step} | ({} if a.action is None else {"action": a.action})
    print(json.dumps({**head, "result": res}, sort_keys=True, default=str))
    return 0


def _build_stage(a) -> int:
    """A stage and its release as a fixed list of run-step steps (card C-04, `biobox.stage`)."""
    from biobox import profile, stage

    libs = tuple(stage.Library(name, version, Path(hmm), cutoff) for name, version, hmm, cutoff in a.library)
    b = stage.Build(root=Path(a.root), stage=a.stage, release=a.release, genomes=Path(a.genomes),
                    objects=Path(a.objects), stages=None if a.stages is None else tuple(_csv(a.stages)),
                    libraries=libs, tmp=Path(a.tmp) if a.tmp else None, threads=a.threads,
                    gene_call_version=a.gene_call_version, scan_jobs=a.scan_jobs, scan_cpus=a.scan_cpus,
                    cluster_parts=a.cluster_parts, ncrna_cm=Path(a.ncrna_cm) if a.ncrna_cm else None,
                    ncrna_set=a.ncrna_set)
    return stage.build_stage(profile.load(a.profile), b)


def main(argv: list[str] | None = None) -> int:
    with nbr.noted_once():  # card W-20: a note about the release (the slow path) once per command, not per store opened
        return _main(argv)


def _main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:2] == ["search", "hmm"]:  # card W-07: its own options, and no recipe or run positional
        return _search_hmm(_hmm_parser().parse_args(argv[2:]))
    command: list[str] = []
    if argv[:1] == ["run-step"] and "--" in argv:  # everything after the first `--` is the step command
        i = argv.index("--")
        argv, command = argv[:i], argv[i + 1:]
    p = argparse.ArgumentParser(prog="bio-box")
    sub = p.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--root", required=True)
    common.add_argument("--gene-call-version", default=GENE_CALL_VERSION)
    pk = sub.add_parser("pack", parents=[common])
    pk.add_argument("action", choices=["build", "verify"])
    pk.add_argument("--release", required=True)
    pk.add_argument("--model", action="append", default=[],
                    help="<library>-<version>=<HMM file>, once for each domain library of the release; optional: "
                         "give none when the release has none (an hmm_file or sequences seed needs no model)")
    pk.add_argument("--pack-dir", help=PACK_DIR_HELP)
    pk.add_argument("--profile", help=PROFILE_DIRS_HELP)
    s = sub.add_parser("search", parents=[common])
    s.add_argument("action", choices=["run", "status", "cancel", "hits"],
                   help="see also `bio-box search hmm --help`: an HMM file to hits and neighbour tables in one command")
    s.add_argument("target", help="recipe path (run), or the run: its directory or its run id in --runs")
    s.add_argument("--release")
    s.add_argument("--wait", action="store_true", help="run: run here until the run closes; its pending record "
                                                         "lets `status` and `cancel` find it meanwhile (card W-13)")
    s.add_argument("--submitted", action="store_true", help=argparse.SUPPRESS)  # the runner `submit` starts
    s.add_argument("--run-id")
    s.add_argument("--work", help="run: parent of the run's work directory (default: <profile tmp>/search/<user>/"
                                  "work with a profile, else $BIOBOX_SEARCH_WORK or /tmp/atd/search)")
    s.add_argument("--runs", help=RUNS_HELP)
    s.add_argument("--pack-dir", help=PACK_DIR_HELP)
    s.add_argument("--cpus", type=int, help=f"worker cpus (default {DEFAULT_CPUS}, or the profile's cap if lower)")
    s.add_argument("--mem-gb", type=float, help="memory cap in GB (default: the profile's mem_scope, else none)")
    s.add_argument("--profile", help="host profile YAML (default: $BIOBOX_PROFILE); caps cpus and memory, "
                                     "and refuses to run uncapped where it requires a cap")
    s.add_argument("--cpu-min", type=float, default=300.0)
    s.add_argument("--wall-s", type=float, default=6 * 3600.0)
    s.add_argument("--round", type=int)
    s.add_argument("--fault", help=argparse.SUPPRESS)  # test hook: oom, timeout, cpu, crash
    s.add_argument("--members", action="store_true")
    s.add_argument("--json", action="store_true")
    s.add_argument("--top_n", type=int, default=export.TOP_N)
    s.add_argument("--format", choices=["tsv", "faa"], help="hits: tsv (one row per gene) or faa (one record per "
                                                           "protein, named by protein uid)")
    s.add_argument("--out", help="run: the run's directory, outside the data root (default: run=<run_id> in "
                                 "--runs); hits --format: write this file (default: stdout)")
    s.add_argument("--library", help="hits --format: domain library of pfam_architecture (default: the only one, "
                                     "else the pfam one)")
    s.add_argument("--threads", type=int, help=f"hits --format: DuckDB threads (default {export.EXPORT_THREADS})")
    e = sub.add_parser("export", parents=[common])
    e.add_argument("what", choices=["proteins"], help="proteins: every distinct protein as FASTA named by protein uid")
    e.add_argument("--release", required=True, help="needs manifests/release=<R>.json (bio-box release manifest)")
    e.add_argument("--out", help="write this file (default: stdout)")
    e.add_argument("--threads", type=int, help=f"DuckDB threads (default {export.EXPORT_THREADS})")
    fd = sub.add_parser("find", parents=[common], help="the genes an accession or a protein uid names (card W-18)")
    fd.add_argument("--release", required=True)
    what = fd.add_mutually_exclusive_group(required=True)
    what.add_argument("--accession", help="a contig accession, or ACCESSION:START-END for the gene with those "
                                          "coordinates (0-based, end exclusive, as `search hits --format tsv` "
                                          "prints), or ACCESSION:START..END (GenBank: 1-based, inclusive)")
    what.add_argument("--protein", help="a protein uid (32 hex digits): every gene of the release that carries it")
    fd.add_argument("--tape-version", help="tapes to read for --accession (default: the gene-call producer's)")
    fd.add_argument("--threads", type=int, help=f"DuckDB threads (default {export.EXPORT_THREADS})")
    rs = sub.add_parser("releases", parents=[common], help="every release of the root: stages, counts, libraries, "
                                                            "lookup tables and pack (card W-23)")
    rs.add_argument("--profile", help="host profile YAML (default: $BIOBOX_PROFILE); with one, the pack of each "
                                      "release in the profile's pack directory is named")
    rs.add_argument("--pack-dir", help="a pack directory to name each release's pack from, in place of the profile's")
    rs.add_argument("--json", action="store_true", help="one JSON object per release")
    sq = sub.add_parser("seq", parents=[common], help="proteins by uid or a contig slice, as FASTA, read by key "
                                                      "(card W-23)")
    sq.add_argument("--release", required=True)
    which = sq.add_mutually_exclusive_group(required=True)
    which.add_argument("--protein", nargs="+", action="extend", help="protein uids (32 hex digits); repeatable")
    which.add_argument("--proteins-file", help="protein uids, one per line (blank lines and # comments skipped)")
    which.add_argument("--contig", help="ACC, ACC:START-END (0-based, end exclusive) or ACC:START..END (GenBank: "
                                        "1-based, inclusive); ACC may be a contig uid. An end past the contig's "
                                        "length crosses the origin of a circular contig")
    sq.add_argument("--strand", type=int, choices=[1, -1], default=1, help="--contig: 1 or -1 (reverse complement)")
    sq.add_argument("--out", help="write this file (default: stdout)")
    sq.add_argument("--threads", type=int, help=f"DuckDB threads (default {export.EXPORT_THREADS})")
    ql = sub.add_parser("sql", parents=[common], help="SELECT statements on a release's catalog, read-only, with the "
                                                      "row groups each read (card W-23)")
    ql.add_argument("--release", required=True)
    ql.add_argument("file", nargs="?", help="a file of SQL statements separated by ;")
    ql.add_argument("-e", "--execute", dest="query", help="the SQL statements, in place of FILE")
    ql.add_argument("--runs", nargs="+", action="extend", help="closed search runs to attach (search_run, search_hit, "
                                                               "...): run directories or runs directories; default "
                                                               "with a profile: its runs directory. Runs of other "
                                                               "releases are left out")
    ql.add_argument("--profile", help="host profile YAML (default: $BIOBOX_PROFILE): its tmp holds the scratch "
                                      "directory, and its runs directory is attached")
    ql.add_argument("--format", choices=["tsv", "json"], default="tsv", help="tsv: a table per statement; json: one "
                                                                             "object per row (one statement)")
    ql.add_argument("--out", help="write the results to this file, outside the data root (default: stdout)")
    ql.add_argument("--threads", type=int, default=4, help="DuckDB threads (default 4)")
    ql.add_argument("--memory-limit", default="4GB", help="DuckDB memory limit; larger work spills into the scratch "
                                                          "directory (default 4GB)")
    n = sub.add_parser("nbr", parents=[common])
    n.add_argument("action", nargs="?", help="summarize: families and clusters per window; flanks: loci grouped by "
                                             "their flanks (card W-25); otherwise the anchors")
    n.add_argument("source", nargs="?", help="summarize, flanks: run_id or protein uids (uid,uid,...)")
    n.add_argument("--anchors", help="run_id or protein uids (uid,uid,...)")
    n.add_argument("--anchors-file", help="protein uids, one per line (blank lines and # comments skipped), instead "
                                          "of --anchors or the uid list")
    n.add_argument("--window", help="10g (genes) or 20kb (bases) on each side; default 10g, flanks: 5g (genes only)")
    n.add_argument("--level", help="flanks: the token id, nr50 (default), nr90 or pfam (the protein's architecture)")
    n.add_argument("--group-by", help="flanks: anchor-nr90 (default), anchor-nr50 or none")
    n.add_argument("--max-loci", type=int, help="flanks: leaves per group in the tree after collapsing identical "
                                                "loci (default 2000; above it, the leaves with the least feature_uid)")
    n.add_argument("--cut", type=float, help="flanks: the tree is cut into context types at this distance "
                                             "(default 0.5)")
    n.add_argument("--distance", help="flanks: the tree's distance, j (default: 1 - Jaccard of the tokens), a "
                                      "(1 - adjacency index) or combined")
    n.add_argument("--export", help="flanks: gff writes one GFF3 per locus and a protein-to-group table")
    n.add_argument("--figure", action="store_true", help="flanks: also groups/<group>/figure.svg and figure.pdf, "
                                                         "the synteny figure (card W-27; needs matplotlib)")
    n.add_argument("--figure-rows", type=int, help="flanks --figure: rows per figure (default 60, "
                                                   "flanks.DEFAULT_FIGURE_ROWS); above it, the leaf with the most loci "
                                                   "of each context type, then the leaves with the most loci")
    n.add_argument("--strand", default="any")
    n.add_argument("--operon-gap")
    n.add_argument("--edge", default="truncate")
    n.add_argument("--circular", default="auto")
    n.add_argument("--json", action="store_true")
    n.add_argument("--top_n", type=int, default=0)
    n.add_argument("--tape-version", help="tapes to read (default: a run's pack's, refused when --pack-dir lacks the "
                                          "pack; protein uids: the gene-call producer's)")
    n.add_argument("--by", default="pfam,nr50", help="summarize: levels, some of pfam,nr50,nr90")
    n.add_argument("--top", type=int, default=nbr.SUMMARY_TOP, help="summarize: rows per level (0: all)")
    n.add_argument("--release", help="summarize, flanks, --format tsv: release of the clusters and contig_dim "
                                     "(default: the run's)")
    n.add_argument("--format", choices=["tsv"], help="summarize: tsv (the default without --json); otherwise one "
                                                     "table per anchor locus")
    n.add_argument("--out", help="--format tsv without summarize: one file per locus in this directory; summarize: "
                                 "write the summary to this file (default: stdout); flanks: a new or empty directory")
    n.add_argument("--library", help="--format tsv without summarize, flanks: domain library of pfam_architecture")
    n.add_argument("--threads", type=int, help=f"--format tsv without summarize, flanks: DuckDB threads (default "
                                               f"{export.EXPORT_THREADS})")
    n.add_argument("--runs", help=RUNS_HELP)
    n.add_argument("--pack-dir", help=PACK_DIR_HELP)
    n.add_argument("--profile", help=PROFILE_DIRS_HELP)
    r = sub.add_parser("run-step", usage="bio-box run-step --profile P --root R --stage S --step NAME "
                                         "[--nice monolithic|fanout] [--threads N] -- CMD [ARGS ...]")
    r.add_argument("--profile", required=True, help="host profile YAML (profiles/wilkcomp02.yaml)")
    r.add_argument("--root", required=True)
    r.add_argument("--stage", required=True)
    r.add_argument("--step", required=True)
    r.add_argument("--nice", choices=["monolithic", "fanout"], default="monolithic")
    r.add_argument("--threads", type=int)
    rl = sub.add_parser("release", parents=[common], help="release-level steps over the release's stages (B-25)")
    rl.add_argument("step", choices=RELEASE_STEPS)
    rl.add_argument("action", nargs="?", choices=CLUSTER_ACTIONS,
                    help="cluster in parts (card B-41): plan, then step (one call per --name, in the plan's "
                         "sequence), then finish; redo makes a step whose scratch was lost again; without it, "
                         "the whole clustering in one process")
    rl.add_argument("--release", required=True)
    rl.add_argument("--name", help="cluster step|redo: the step (cluster.plan, cluster.<level>.<kk>.<module>, its "
                                   ".cut or .part-<p>, ...)")
    rl.add_argument("--parts", type=int, help="cluster plan: cut every per-query call into this many parts (forced)")
    rl.add_argument("--part-cpu-min", type=float, help="cluster plan: the CPU-minutes a part may be estimated at "
                                                       "(build-stage gives 80%% of cpu_min_per_process)")
    rl.add_argument("--one-process-cpu-min-max", type=float,
                    help="cluster plan: the cap under which a call that no number of parts fits runs as one process "
                         "(DEC-29, card B-45; build-stage gives the profile's cluster_cpu_min_max); without it such a "
                         "call is refused")
    rl.add_argument("--stages", help="S0,S1: the release's stages in order; defines the release once")
    rl.add_argument("--tmp", help="scratch directory (system temp by default; never under the data root)")
    rl.add_argument("--threads", type=int, default=1, help="cluster: MMseqs2 threads")
    rl.add_argument("--producers", help="tokens, postings: the tape version's producers (default: the gene caller)")
    rl.add_argument("--libraries", help="tokens: domain libraries in precedence order (default: those with domains in the release's stages, by name)")
    rl.add_argument("--tape-version", help="catalog: the tape version the pair macros serve (default: the one of "
                                           "the release's postings)")
    rl.add_argument("--out", help=OUT_HELP)
    rl.add_argument("--previous", help="verify: an earlier release; every file of its manifest outside its own "
                                       "release= partition must be in this release's with the same sha256 (R-09)")
    st = sub.add_parser("stage", parents=[common], help="per-stage builder steps (card C-08)")
    st.add_argument("step", choices=STAGE_STEPS)
    st.add_argument("action", nargs="?", choices=tuple(dict.fromkeys(DOMAIN_ACTIONS + INGEST_ACTIONS + NCRNA_ACTIONS)),
                    help="domains (card B-24) or ncrna (card E-07) in shards: plan, then shard (one call per --shard "
                         "K), then merge; "
                         "ingest in batch processes (card B-30): pass1 (one call per --shard K), plan, batch (one "
                         "call per --batch J), then merge; without it, the whole step in one process")
    st.add_argument("--stage", required=True)
    st.add_argument("--out", help=OUT_HELP)
    st.add_argument("--tmp", help="scratch directory (system temp by default; under neither the root nor --out)")
    st.add_argument("--max-file-bytes", type=int, help="cut each table into parts of at most this many estimated bytes")
    st.add_argument("--genomes", help="inputs: JSON list of genomes (GenomeInput fields, contig uids as hex)")
    st.add_argument("--objects", help="ingest: JSON list of source objects (SourceObject fields)")
    st.add_argument("--shard-bases", type=int, help="ingest: new bases per sequence shard")
    st.add_argument("--line-width", type=int, help="ingest: FASTA line width of the shards")
    st.add_argument("--prefix-width", type=int, choices=[0, 1, 2], help="genecall: protein_seq bucket width")
    st.add_argument("--workers", type=int, default=1, help="genecall: genomes called at once")
    st.add_argument("--seq-stages", help="genecall, ncrna: S0,S1: the stages whose sequences the genomes name "
                                         "(default: the stage itself)")
    st.add_argument("--producers", help="igr, tapes: the tape version's producers (default: the gene caller)")
    st.add_argument("--max-slots", type=int, help="tapes: slots per row group")
    st.add_argument("--library", help="domains: library name (pfam)")
    st.add_argument("--library-version", help="domains: library version (38.2)")
    st.add_argument("--hmm", help="domains: the library's HMM file")
    st.add_argument("--cutoff", help="domains: gathering, trusted, noise or none")
    st.add_argument("--cpus", type=int, default=1, help="domains: hmmsearch threads; ncrna: cmsearch threads (--cpu)")
    st.add_argument("--cm", help="ncrna: the Rfam covariance-model library (.cm or .cm.gz; card E-01)")
    st.add_argument("--ncrna-set", default="prok-2", help="ncrna: the model set, a frozen set's name (prok-2, the "
                                                          "default, DEC-30; prok-1) or a set file (card E-07)")
    st.add_argument("--shard-residues", type=int, help="domains plan: residues per shard at most")
    st.add_argument("--scan-bases", type=int, help="ncrna plan: genome bases per shard at most")
    st.add_argument("--shard", type=int, help="domains or ncrna shard: the shard to scan; ingest pass1: the pass-1 "
                                             "shard (0-based)")
    st.add_argument("--pass1-bytes", type=int, help="ingest pass1: bytes of source file per pass-1 shard at most")
    st.add_argument("--batch", type=int, help="ingest batch: the batch to write (0-based)")
    st.add_argument("--runs-tmp", help="ingest pass1|plan|batch|merge: the scratch directory the steps share for "
                                       "the pass-1 runs and the decision table (default: --tmp)")
    bs = sub.add_parser("build-stage", parents=[common], help="a stage and its release, step by step under "
                                                              "run-step (card C-04)")
    bs.add_argument("--profile", required=True, help="host profile YAML (profiles/wilkcomp02.yaml)")
    bs.add_argument("--stage", required=True)
    bs.add_argument("--release", required=True)
    bs.add_argument("--stages", help="S0,S1: the release's stages in order; defines the release on first use, "
                                     "else they come from its release_def")
    bs.add_argument("--genomes", required=True, help="inputs: JSON list of genomes (as bio-box stage inputs)")
    bs.add_argument("--objects", required=True, help="ingest: JSON list of source objects (as bio-box stage ingest)")
    bs.add_argument("--library", nargs=4, action="append", default=[], metavar=("NAME", "VERSION", "HMM", "CUTOFF"),
                    help="a domain library to scan (repeatable; one step each, in the order given)")
    bs.add_argument("--tmp", help="scratch directory of the steps (default: <profile tmp>/scratch)")
    bs.add_argument("--threads", type=int, help="threads of every step but the scan shards (default: the smaller "
                                                "of the profile's day and night caps)")
    bs.add_argument("--scan-jobs", type=int, help="scan shards run at once, each its own run-step (card C-11; "
                                                  "default: the threads run-step allows now // --scan-cpus)")
    bs.add_argument("--scan-cpus", type=int, help="threads of each scan shard (default 4, or the allowance if "
                                                  "smaller); --scan-jobs × --scan-cpus may not exceed the threads "
                                                  "run-step allows now")
    bs.add_argument("--cluster-parts", type=int, help="cluster in parts (card B-41), every per-query call cut into "
                                                      "this many (default: in parts only when clustering in one "
                                                      "process is estimated over 80%% of the budget, cut by it)")
    bs.add_argument("--ncrna-cm", help="the Rfam library (Rfam.cm.gz): run the ncRNA producer in the features slot and "
                                       "build IGRs, tapes, tokens and postings of both producers (card E-07)")
    bs.add_argument("--ncrna-set", default="prok-2", help="the ncRNA model set: a frozen set's name (prok-2, the "
                                                          "default, DEC-30; prok-1) or a set file")
    a = p.parse_args(argv)
    if a.cmd == "search" and a.action == "run" and not a.release:
        p.error("search run needs --release")
    if a.cmd == "run-step":
        if not command:
            p.error("run-step needs a command after --")
        a.command = command
    return {"pack": _pack, "search": _search, "nbr": _nbr, "find": _find, "export": _export, "run-step": _run_step,
            "release": _release, "stage": _stage, "build-stage": _build_stage, "releases": _releases,
            "seq": _seq, "sql": _sql}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
