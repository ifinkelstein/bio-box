"""`bio-box` command line: search pack, search runs, neighborhoods (SEARCH_PLAN §3.3).

    bio-box pack build --root R --release rS0.1 [--model pfam-38.2=/path/Pfam-A.hmm ...] [--pack-dir P]
    bio-box pack verify --root R --release rS0.1 [--pack-dir P]
    bio-box search run recipe.yaml --root R --release rS0.1 [--wait] [--out DIR | --runs DIR] [--pack-dir P]
    bio-box search status|cancel <run> --root R [--runs DIR]
    bio-box search hits <run> --root R [--runs DIR] [--round N] [--members] [--json] [--top_n 50]
    bio-box search hits <run> --root R [--runs DIR] --format tsv|faa [--round N] [--library L] [--out FILE]
    bio-box nbr --root R --anchors <run|uid,uid,...> --window 10g --strand any --edge truncate [--json]
                [--tape-version T] [--runs DIR]
    bio-box nbr <run|uid,uid,...> --root R --window 5g --format tsv [--release R] [--out DIR] [--runs DIR]
    bio-box nbr summarize <run|uid,uid,...> --root R --window 5g --by pfam,nr50 [--top 50] [--release R]
                [--json | --format tsv] [--runs DIR]
    bio-box export proteins --root R --release R [--out FILE]
    bio-box run-step --profile P --root R --stage S --step NAME [--nice monolithic|fanout] [--threads N] -- CMD ...
    bio-box release define --root R --release R1 --stages S0,S1 [--out DIR]
    bio-box release cluster|tokens|postings|dim|manifest|catalog --root R --release R1 [--stages S0,S1]
                [--tmp D] [--threads N] [--producers P,Q] [--libraries L,M] [--tape-version T] [--out DIR]
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
    bio-box stage ncrna --root R --stage S --cm Rfam.cm.gz --cpus N [--tmp D] [--seq-stages S0,S1] [--out DIR]
    bio-box build-stage --profile P --root R --stage S --release R1 [--stages S0,S1] --genomes genomes.json
                --objects objects.json [--library NAME VERSION HMM CUTOFF ...] [--tmp D] [--threads N]
                [--scan-jobs J] [--scan-cpus C]

`build-stage` (card C-04, `biobox.stage`) runs every `stage` step and then every `release` step of
one stage and its release, each as a `run-step`; a second call skips what is done.

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
missing exits 2 with the command that builds it.

`search run --wait` runs the search in this process and writes the same pending record as a submitted
run, so `search status` and `search cancel` work on it from another shell (card W-13). Without `--wait`
the run is submitted: a detached runner, and the command returns at once. Inside a `scripts/pod_run.sh`
unit (BIOBOX_POD_UNIT, which the script sets) that is refused (exit 2): the unit is a transient service
started with `systemd-run --collect`, which ends when its main process exits and stops every process
left in its cgroup, the detached runner included; a unit stays alive only as long as its main process
runs, and with `--wait` the main process is the runner.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq

from biobox import catalog
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.search import export, nbr, pack, runner
from biobox.search.recipe import RecipeError


def _print(lines) -> None:
    for line in lines:
        sys.stdout.write(line + "\n")


def _emit(lines, out: str | None) -> None:
    """To stdout, or to the file `out` (written under a temporary name and renamed into place)."""
    if out is None:
        _print(lines)
    else:
        export.write_lines(lines, Path(out))


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


def _run_tape_version(d: Path, packs: Path, gene_call_version: str) -> str:
    """The tape version a closed run read: its pack's (card B-16), the pack found by the run's pack_id in
    the pack directory `packs`; the gene-call producer's when the run or its pack is not on disk."""
    path = Path(d) / "run.parquet"
    if path.exists():
        row = pq.read_table(path, columns=["ran_against", "pack_id"]).to_pylist()[0]
        pdir = pack.pack_dir(packs, row["ran_against"], row["pack_id"])
        if (pdir / pack.MANIFEST).exists():
            return pack.manifest_tape_version(pack.read_manifest(pdir)[0])
    return nbr.resolve_tape_version(gene_call_version)


def _pack(a) -> int:
    """`pack build|verify` (card W-03): the release's pack in the pack directory, never under the root."""
    packs = _packs(a)
    if a.action == "build":
        models = dict(m.split("=", 1) for m in a.model)
        try:
            man = pack.build(Path(a.root), a.release, a.gene_call_version, {k: Path(v) for k, v in models.items()},
                             packs=packs)
        except (ValueError, FileExistsError) as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        recs = (export.record("pack_target", {"release_id": a.release, "pack_id": man["pack_id"], "target": t,
                                              "shards": sum(1 for s in man["shards"] if s["target"] == t), **info})
                for t, info in sorted(man["targets"].items()))
        _print(export.ndjson(recs))
        return 0
    problems = pack.verify(Path(a.root), a.release, packs=packs)
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
        cap = runstep.choose_threads(prof, datetime.now(), runstep._load1())
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
            print(f"refused: no pending run {a.target} (looked for {runner.pending_path(d)})", file=sys.stderr)
            return 2
        return 0
    if not (d / "run.parquet").exists():
        print(f"refused: no closed run {a.target} (looked for {d})", file=sys.stderr)
        return 2
    hits = runner.read_hits(d, a.round)
    tv = _run_tape_version(d, _packs(a), a.gene_call_version)
    release = _run_release(d)
    manifest = export.release_manifest(root, release)
    if a.format:  # card W-05 (files: the manifest's, else the release_def stages: card W-12): one row per gene (tsv) or one record per protein (faa)
        if a.json:
            print("refused: --json and --format are exclusive", file=sys.stderr)
            return 2
        make = export.hits_tsv if a.format == "tsv" else export.hits_faa
        try:
            _emit(make(root, hits, release=release, tape_version=tv, gene_call_version=a.gene_call_version,
                       manifest=manifest, library=a.library, threads=a.threads), a.out)
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
    return 0


def _anchor_run(a, source: str) -> Path | None:
    """The run directory `source` names (`_run_path`), when it holds a closed run; else None (protein uids)."""
    d = _run_path(a, source)
    return d if (d / "run.parquet").exists() else None


def _anchors(a, source: str, tape_version: str | None, gene_call_version: str):
    """Anchor proteins from a recorded run's protein hits (a run directory, or a run id in the runs
    directory) or a comma-separated list of protein uids: (uids, tape version, the run's release or None),
    or a refusal message."""
    d = _anchor_run(a, source)
    if d is not None:
        uids = {h["target_uid"] for h in runner.read_hits(d) if h["target_kind"] == "protein"}
        run_tv = _run_tape_version(d, _packs(a), gene_call_version)
        if tape_version is not None and tape_version != run_tv:
            return f"refused: tape version {tape_version} is not the tape version {run_tv} that run {source} read"
        return uids, run_tv, _run_release(d)
    try:
        return {bytes.fromhex(x) for x in source.split(",") if x}, tape_version, None
    except ValueError:
        return f"refused: {source} is neither a closed run (in {_runs(a)} or as a directory) nor protein uids"


def _nbr(a) -> int:
    root = Path(a.root)
    if a.action not in (None, "summarize"):  # `nbr <run_id|uids>`: the first positional is the anchors
        if a.source is not None:
            print(f"refused: nbr: unknown action {a.action!r} (expected summarize)", file=sys.stderr)
            return 2
        a.action, a.source = None, a.action
    if a.action == "summarize" and not (a.source or a.anchors):
        print("refused: nbr summarize needs a run id or a list of protein uids", file=sys.stderr)
        return 2
    if a.action is None and not (a.anchors or a.source):
        print("refused: nbr needs --anchors", file=sys.stderr)
        return 2
    source = a.source or a.anchors
    got = _anchors(a, source, a.tape_version, a.gene_call_version)
    if isinstance(got, str):
        print(got, file=sys.stderr)
        return 2
    uids, tv, run_release = got
    policy = nbr.Policy.parse(a.window, a.strand, a.operon_gap, a.edge, a.circular)
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
                for name, lines in tables:
                    export.write_lines(lines, Path(a.out) / name)
        except ValueError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        return 0
    if a.action == "summarize":
        try:
            s = nbr.summarize(root, uids, policy, [x.strip() for x in a.by.split(",") if x.strip()], release=release,
                              gene_call_version=a.gene_call_version, tape_version=tv, top_n=a.top,
                              manifest=files.manifest)
        except ValueError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        _print(export.ndjson(export.nbr_summary_records(s)) if a.json else export.nbr_summary_tsv(s))
        return 0
    rows, counts = nbr.fetch(root, uids, policy, a.gene_call_version, tape_version=tv, files=files)
    if a.json:
        _print(export.ndjson(export.nbr_records(rows, a.top_n or None)))
    else:
        print(json.dumps(counts))
    return 0


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
                               threads=a.threads), a.out)
    return 0


def _run_step(a) -> int:
    from biobox import profile, runstep

    return runstep.run_step(profile.load(a.profile), Path(a.root), a.stage, a.step, a.command,
                            nice=a.nice, threads=a.threads)


RELEASE_STEPS = ("define", "cluster", "tokens", "postings", "dim", "manifest", "catalog", "verify")
STAGE_STEPS = ("inputs", "ingest", "genecall", "igr", "tapes", "domains", "ncrna")
DOMAIN_ACTIONS = ("plan", "shard", "merge")  # `stage domains plan|shard|merge` (card B-24)
INGEST_ACTIONS = ("pass1", "plan", "batch", "merge")  # `stage ingest pass1|plan|batch|merge` (card B-30)
STEP_ACTIONS = {"domains": DOMAIN_ACTIONS, "ingest": INGEST_ACTIONS}
PACK_DIR_HELP = ("the pack directory, outside the data root (card W-03); default: <profile tmp>/packs with a host "
                 "profile, else packs beside $BIOBOX_SEARCH_WORK or /tmp/atd/search")
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
    when there is no problem, else 1 (a missing manifest is a problem)."""
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
    failed = bool(problems)
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
    print(json.dumps({"release": a.release, "stages": stages, "step": a.step, "result": res}, sort_keys=True,
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
    if a.step == "ncrna":  # card E-01
        from biobox.l1 import ncrna

        if not a.cm:
            raise ValueError("stage ncrna needs --cm")
        return ncrna.run(root, a.stage, Path(a.cm), a.cpus, tmp=tmp, out=out, seq_stages=_csv(a.seq_stages), **mfb)
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
                    gene_call_version=a.gene_call_version, scan_jobs=a.scan_jobs, scan_cpus=a.scan_cpus)
    return stage.build_stage(profile.load(a.profile), b)


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
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
    pk.add_argument("--model", action="append", default=[], help="library=path to an HMM file")
    pk.add_argument("--pack-dir", help=PACK_DIR_HELP)
    pk.add_argument("--profile", help=PROFILE_DIRS_HELP)
    s = sub.add_parser("search", parents=[common])
    s.add_argument("action", choices=["run", "status", "cancel", "hits"])
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
    n = sub.add_parser("nbr", parents=[common])
    n.add_argument("action", nargs="?", help="summarize: families and clusters per window; otherwise the anchors")
    n.add_argument("source", nargs="?", help="summarize: run_id or protein uids (uid,uid,...)")
    n.add_argument("--anchors", help="run_id or protein uids (uid,uid,...)")
    n.add_argument("--window", default="10g")
    n.add_argument("--strand", default="any")
    n.add_argument("--operon-gap")
    n.add_argument("--edge", default="truncate")
    n.add_argument("--circular", default="auto")
    n.add_argument("--json", action="store_true")
    n.add_argument("--top_n", type=int, default=0)
    n.add_argument("--tape-version", help="tapes to read (default: the run's pack's, else the gene-call producer's)")
    n.add_argument("--by", default="pfam,nr50", help="summarize: levels, some of pfam,nr50,nr90")
    n.add_argument("--top", type=int, default=nbr.SUMMARY_TOP, help="summarize: rows per level (0: all)")
    n.add_argument("--release", help="summarize, --format tsv: release of the clusters and contig_dim (default: "
                                     "the run's)")
    n.add_argument("--format", choices=["tsv"], help="summarize: tsv (the default without --json); otherwise one "
                                                     "table per anchor locus")
    n.add_argument("--out", help="--format tsv without summarize: one file per locus in this directory")
    n.add_argument("--library", help="--format tsv without summarize: domain library of pfam_architecture")
    n.add_argument("--threads", type=int, help=f"--format tsv without summarize: DuckDB threads (default "
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
    rl.add_argument("--release", required=True)
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
    st.add_argument("action", nargs="?", choices=tuple(dict.fromkeys(DOMAIN_ACTIONS + INGEST_ACTIONS)),
                    help="domains in shards (card B-24): plan, then shard (one call per --shard K), then merge; "
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
    st.add_argument("--seq-stages", help="genecall: S0,S1: the stages whose sequences the genomes name "
                                         "(default: the stage itself)")
    st.add_argument("--producers", help="igr, tapes: the tape version's producers (default: the gene caller)")
    st.add_argument("--max-slots", type=int, help="tapes: slots per row group")
    st.add_argument("--library", help="domains: library name (pfam)")
    st.add_argument("--library-version", help="domains: library version (38.2)")
    st.add_argument("--hmm", help="domains: the library's HMM file")
    st.add_argument("--cutoff", help="domains: gathering, trusted, noise or none")
    st.add_argument("--cpus", type=int, default=1, help="domains: hmmsearch threads; ncrna: cmsearch threads (--cpu)")
    st.add_argument("--cm", help="ncrna: the Rfam covariance-model library (.cm or .cm.gz; card E-01)")
    st.add_argument("--shard-residues", type=int, help="domains plan: residues per shard at most")
    st.add_argument("--shard", type=int, help="domains shard: the shard to scan; ingest pass1: the pass-1 shard "
                                             "(0-based)")
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
    a = p.parse_args(argv)
    if a.cmd == "search" and a.action == "run" and not a.release:
        p.error("search run needs --release")
    if a.cmd == "run-step":
        if not command:
            p.error("run-step needs a command after --")
        a.command = command
    return {"pack": _pack, "search": _search, "nbr": _nbr, "export": _export, "run-step": _run_step,
            "release": _release, "stage": _stage, "build-stage": _build_stage}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
