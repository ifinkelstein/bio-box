"""`bio-box` command line: search pack, search runs, neighborhoods (SEARCH_PLAN §3.3).

    bio-box pack build --root R --release rS0.1 [--model pfam-38.2=/path/Pfam-A.hmm ...]
    bio-box pack verify --root R --release rS0.1
    bio-box search run recipe.yaml --root R --release rS0.1 [--wait]
    bio-box search status|cancel <run_id> --root R
    bio-box search hits <run_id> --root R [--round N] [--members] [--json] [--top_n 50]
    bio-box search hits <run_id> --root R --format tsv|faa [--round N] [--library L] [--out FILE]
    bio-box nbr --root R --anchors <run_id|uid,uid,...> --window 10g --strand any --edge truncate [--json]
                [--tape-version T]
    bio-box nbr <run_id|uid,uid,...> --root R --window 5g --format tsv [--release R] [--out DIR]
    bio-box nbr summarize <run_id|uid,uid,...> --root R --window 5g --by pfam,nr50 [--top 50] [--release R]
                [--json | --format tsv]
    bio-box export proteins --root R --release R [--out FILE]
    bio-box run-step --profile P --root R --stage S --step NAME [--nice monolithic|fanout] [--threads N] -- CMD ...
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq

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


def _run_release(root: Path, run_id: str) -> str | None:
    """The release a recorded run ran against; None when the run is not on disk."""
    path = runner.run_dir(root, run_id) / "run.parquet"
    if not path.exists():
        return None
    return pq.read_table(path, columns=["ran_against"]).column("ran_against")[0].as_py()


def _run_tape_version(root: Path, run_id: str, gene_call_version: str) -> str:
    """The tape version a closed run read: its pack's (card B-16); the gene-call producer's when the
    run or its pack is not on disk."""
    path = runner.run_dir(root, run_id) / "run.parquet"
    if path.exists():
        release = pq.read_table(path, columns=["ran_against"]).column("ran_against")[0].as_py()
        if (pack.pack_dir(root, release) / "pack_manifest.json").exists():
            return pack.manifest_tape_version(pack.read_manifest(root, release)[0])
    return nbr.resolve_tape_version(gene_call_version)


def _pack(a) -> int:
    if a.action == "build":
        models = dict(m.split("=", 1) for m in a.model)
        man = pack.build(Path(a.root), a.release, a.gene_call_version, {k: Path(v) for k, v in models.items()})
        recs = (export.record("pack_target", {"release_id": a.release, "pack_id": man["pack_id"], "target": t,
                                              "shards": sum(1 for s in man["shards"] if s["target"] == t), **info})
                for t, info in sorted(man["targets"].items()))
        _print(export.ndjson(recs))
        return 0
    problems = pack.verify(Path(a.root), a.release)
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


def _search(a) -> int:
    root = Path(a.root)
    if a.action == "run":
        try:
            limits = _limits(a)
        except (ValueError, OSError) as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        try:
            if not a.wait:
                print(runner.submit(root, Path(a.target), a.release, limits, Path(a.work) if a.work else None,
                                    a.fault))
                return 0
            row = runner.run(root, Path(a.target), a.release, limits=limits,
                             work_base=Path(a.work) if a.work else None, argv=sys.argv,
                             run_id=bytes.fromhex(a.run_id) if a.run_id else None, fault=a.fault)
        except RecipeError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        _print(export.ndjson([export.run_record(row)]) if a.json else [row["run_id"].hex() + "\t" + row["status"]])
        return 0 if row["status"] == "ok" else 1
    if a.action == "status":
        print(json.dumps(runner.status(root, a.target)))
        return 0
    if a.action == "cancel":
        runner.cancel(root, a.target)
        return 0
    hits = runner.read_hits(root, a.target, a.round)
    tv = _run_tape_version(root, a.target, a.gene_call_version)
    release = _run_release(root, a.target)
    manifest = export.release_manifest(root, release)
    if a.format:  # card W-05: one row per gene (tsv) or one record per protein (faa)
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
    for o in nbr.occurrences(root, {h["target_uid"] for h in hits if h["target_kind"] == "protein"},
                             a.gene_call_version, tape_version=tv):
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


def _anchors(root: Path, source: str, tape_version: str | None, gene_call_version: str):
    """Anchor proteins from a recorded run's protein hits or a comma-separated list of protein uids:
    (uids, tape version, the run's release or None), or a refusal message."""
    if len(source) == 32 and (runner.run_dir(root, source) / "run.parquet").exists():
        uids = {h["target_uid"] for h in runner.read_hits(root, source) if h["target_kind"] == "protein"}
        run_tv = _run_tape_version(root, source, gene_call_version)
        if tape_version is not None and tape_version != run_tv:
            return f"refused: tape version {tape_version} is not the tape version {run_tv} that run {source} read"
        release = pq.read_table(runner.run_dir(root, source) / "run.parquet",
                                columns=["ran_against"]).column("ran_against")[0].as_py()
        return uids, run_tv, release
    return {bytes.fromhex(x) for x in source.split(",") if x}, tape_version, None


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
    got = _anchors(root, source, a.tape_version, a.gene_call_version)
    if isinstance(got, str):
        print(got, file=sys.stderr)
        return 2
    uids, tv, run_release = got
    policy = nbr.Policy.parse(a.window, a.strand, a.operon_gap, a.edge, a.circular)
    if a.action == "summarize" or a.format == "tsv":
        release = a.release or run_release
        if release is None:
            print("refused: anchors given as protein uids need --release (the release whose clusters and "
                  "contig_dim to read)", file=sys.stderr)
            return 2
    if a.action is None and a.format == "tsv":  # card W-05: one table per anchor locus
        if a.json:
            print("refused: --json and --format are exclusive", file=sys.stderr)
            return 2
        rows, _ = nbr.fetch(root, uids, policy, a.gene_call_version, tape_version=tv)
        bits = None
        if run_release is not None:
            bits = {}
            for h in runner.read_hits(root, source):
                if h["target_kind"] == "protein" and h.get("seq_bitscore") is not None:
                    bits[h["target_uid"]] = max(bits.get(h["target_uid"], h["seq_bitscore"]), h["seq_bitscore"])
        try:
            tables = export.locus_tables(root, rows, release=release, gene_call_version=a.gene_call_version,
                                         hit_bits=bits, manifest=export.release_manifest(root, release),
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
                              gene_call_version=a.gene_call_version, tape_version=tv, top_n=a.top)
        except ValueError as e:
            print(f"refused: {e}", file=sys.stderr)
            return 2
        _print(export.ndjson(export.nbr_summary_records(s)) if a.json else export.nbr_summary_tsv(s))
        return 0
    rows, counts = nbr.fetch(root, uids, policy, a.gene_call_version, tape_version=tv)
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
        print(f"note: no manifests/release={a.release}.json; reading every protein_seq file of gene_call_version "
              f"{a.gene_call_version}", file=sys.stderr)
    _emit(export.protein_fasta(root, release=a.release, gene_call_version=a.gene_call_version, manifest=manifest,
                               threads=a.threads), a.out)
    return 0


def _run_step(a) -> int:
    from biobox import profile, runstep

    return runstep.run_step(profile.load(a.profile), Path(a.root), a.stage, a.step, a.command,
                            nice=a.nice, threads=a.threads)


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
    s = sub.add_parser("search", parents=[common])
    s.add_argument("action", choices=["run", "status", "cancel", "hits"])
    s.add_argument("target", help="recipe path (run) or run_id")
    s.add_argument("--release")
    s.add_argument("--wait", action="store_true")
    s.add_argument("--run-id")
    s.add_argument("--work")
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
    s.add_argument("--out", help="hits --format: write this file (default: stdout)")
    s.add_argument("--library", help="hits --format: domain library of pfam_architecture (default: the only one, "
                                     "else the pfam one)")
    s.add_argument("--threads", type=int, help=f"hits --format: DuckDB threads (default {export.EXPORT_THREADS})")
    e = sub.add_parser("export", parents=[common])
    e.add_argument("what", choices=["proteins"], help="proteins: every distinct protein as FASTA named by protein uid")
    e.add_argument("--release", required=True, help="files from manifests/release=<R>.json when present")
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
    r = sub.add_parser("run-step", usage="bio-box run-step --profile P --root R --stage S --step NAME "
                                         "[--nice monolithic|fanout] [--threads N] -- CMD [ARGS ...]")
    r.add_argument("--profile", required=True, help="host profile YAML (profiles/wilkcomp02.yaml)")
    r.add_argument("--root", required=True)
    r.add_argument("--stage", required=True)
    r.add_argument("--step", required=True)
    r.add_argument("--nice", choices=["monolithic", "fanout"], default="monolithic")
    r.add_argument("--threads", type=int)
    a = p.parse_args(argv)
    if a.cmd == "search" and a.action == "run" and not a.release:
        p.error("search run needs --release")
    if a.cmd == "run-step":
        if not command:
            p.error("run-step needs a command after --")
        a.command = command
    return {"pack": _pack, "search": _search, "nbr": _nbr, "export": _export, "run-step": _run_step}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
