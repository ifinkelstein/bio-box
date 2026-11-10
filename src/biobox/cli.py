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
    bio-box release define --root R --release R1 --stages S0,S1 [--out DIR]
    bio-box release cluster|tokens|postings|dim|manifest|catalog --root R --release R1 [--stages S0,S1]
                [--tmp D] [--threads N] [--producers P,Q] [--libraries L,M] [--tape-version T] [--out DIR]
    bio-box release verify --root R --release R2 [--previous R1]
    bio-box stage inputs --root R --stage S --genomes genomes.json [--out DIR]
    bio-box stage ingest --root R --stage S --objects objects.json [--shard-bases N] [--line-width N] [--tmp D]
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
            for h in runner.read_hits(root, source):
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
    if a.action is not None and a.step != "domains":
        raise ValueError(f"stage {a.step}: {a.action} is a step of stage domains only")
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
        return ingest(_source_objects(a.objects), root, a.stage, tmp=tmp, out=out, **kw, **mfb)
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
    st.add_argument("action", nargs="?", choices=DOMAIN_ACTIONS,
                    help="domains in shards (card B-24): plan, then shard (one call per --shard K), then merge; "
                         "without it, the whole scan in one process")
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
    st.add_argument("--shard", type=int, help="domains shard: the shard to scan (0-based)")
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
