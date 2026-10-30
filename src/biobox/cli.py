"""`bio-box` command line: search pack, search runs, neighborhoods (SEARCH_PLAN §3.3).

    bio-box pack build --root R --release rS0.1 [--model pfam-38.2=/path/Pfam-A.hmm ...]
    bio-box pack verify --root R --release rS0.1
    bio-box search run recipe.yaml --root R --release rS0.1 [--wait]
    bio-box search status|cancel <run_id> --root R
    bio-box search hits <run_id> --root R [--round N] [--members] [--json] [--top_n 50] [--format tsv]
    bio-box nbr --root R --anchors <run_id|uid,uid,...> --window 10g --strand any --edge truncate [--json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.search import export, nbr, pack, runner
from biobox.search.recipe import RecipeError


def _print(lines) -> None:
    for line in lines:
        sys.stdout.write(line + "\n")


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


def _search(a) -> int:
    root = Path(a.root)
    if a.action == "run":
        limits = runner.Limits(cpus=a.cpus, cpu_min=a.cpu_min, wall_s=a.wall_s)
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
    occ = {}
    for o in nbr.occurrences(root, {h["target_uid"] for h in hits if h["target_kind"] == "protein"},
                             a.gene_call_version):
        occ.setdefault(o["protein_uid"], []).append(o)
    recs = list(export.hit_records(hits, occ, a.top_n or None, a.members))
    if a.json:
        _print(export.ndjson(recs))
    else:
        cols = ["target_kind", "target_uid", "seq_bitscore", "bitscore", "evalue_zref", "env_start", "env_end",
                "n_occurrences", "curation_label"]
        _print(["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in recs])
    return 0


def _nbr(a) -> int:
    root = Path(a.root)
    if len(a.anchors) == 32 and (runner.run_dir(root, a.anchors) / "run.parquet").exists():
        uids = {h["target_uid"] for h in runner.read_hits(root, a.anchors) if h["target_kind"] == "protein"}
    else:
        uids = {bytes.fromhex(x) for x in a.anchors.split(",") if x}
    policy = nbr.Policy.parse(a.window, a.strand, a.operon_gap, a.edge, a.circular)
    rows, counts = nbr.fetch(root, uids, policy, a.gene_call_version)
    if a.json:
        _print(export.ndjson(export.nbr_records(rows, a.top_n or None)))
    else:
        print(json.dumps(counts))
    return 0


def main(argv: list[str] | None = None) -> int:
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
    s.add_argument("--cpus", type=int, default=8)
    s.add_argument("--cpu-min", type=float, default=300.0)
    s.add_argument("--wall-s", type=float, default=6 * 3600.0)
    s.add_argument("--round", type=int)
    s.add_argument("--fault", help=argparse.SUPPRESS)  # test hook: oom, timeout, cpu, crash
    s.add_argument("--members", action="store_true")
    s.add_argument("--json", action="store_true")
    s.add_argument("--top_n", type=int, default=export.TOP_N)
    n = sub.add_parser("nbr", parents=[common])
    n.add_argument("--anchors", required=True)
    n.add_argument("--window", default="10g")
    n.add_argument("--strand", default="any")
    n.add_argument("--operon-gap")
    n.add_argument("--edge", default="truncate")
    n.add_argument("--circular", default="auto")
    n.add_argument("--json", action="store_true")
    n.add_argument("--top_n", type=int, default=0)
    a = p.parse_args(argv)
    if a.cmd == "search" and a.action == "run" and not a.release:
        p.error("search run needs --release")
    return {"pack": _pack, "search": _search, "nbr": _nbr}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
